"""Multi-station relay-plan construction.

Given one satellite and several ground stations, the existing
:func:`app.passes.find_passes` proves per-station visibility windows.
This module turns those proven windows into a *relay plan*: a single
longest continuous stretch of the search interval over which the
satellite is provably visible from at least one of the stations and the
stations can hand off to each other.

Only proven visibility is ever used. A handoff from station A to station
B is legal only where *both* windows cover a shared interval of positive
duration; two windows merely touching at one instant cannot hand off.
Time that no station proves visible but at least one station could not
resolve is reported as *unknown* rather than claimed invisible.

Half-open semantics are inherited from :mod:`app.passes`: every window is
``[start, end)``. Positive overlap therefore means
``max(start_a, start_b) < min(end_a, end_b)``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional

from .passes import PassResult


@dataclass
class ProvenWindow:
    """One half-open proven-visible window [start, end) of one station."""
    station: str
    start: datetime
    end: datetime
    truncated: bool


@dataclass
class Assignment:
    """Station ``station`` is responsible for the half-open interval
    ``[start, end)`` of the relay plan."""
    station: str
    start: datetime
    end: datetime
    duration_s: float
    truncated: bool


@dataclass
class Handoff:
    """Responsibility switches from ``from_station`` to ``to_station`` at
    ``time`` (both stations are proven visible on an interval of positive
    duration around the switch)."""
    time: datetime
    from_station: str
    to_station: str


@dataclass
class UnknownInterval:
    start: datetime
    end: datetime
    reasons: list[str]


@dataclass
class RelayPlan:
    start: datetime
    end: datetime
    duration_s: float
    switches: int
    assignments: list[Assignment] = field(default_factory=list)
    handoffs: list[Handoff] = field(default_factory=list)
    stations: tuple[str, ...] = ()


@dataclass
class StationSummary:
    id: str
    status: str
    windows: list[ProvenWindow] = field(default_factory=list)
    uncovered: list = field(default_factory=list)  # list[UncoveredInterval]
    truncated: bool = False


@dataclass
class RelayResult:
    status: str                       # complete | partial | failed
    plan: Optional[RelayPlan]
    stations: list[StationSummary]
    unknown_intervals: list[UnknownInterval]


# --------------------------------------------------------------- intervals
def _merge(intervals):
    """Union of (a, b) intervals; intervals merely touching at an instant
    stay separate (a half-open relay needs positive-length coverage)."""
    if not intervals:
        return []
    intervals = sorted(intervals)
    out = [list(intervals[0])]
    for a, b in intervals[1:]:
        if a < out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]


def _subtract(base, holes):
    """Set difference of interval lists."""
    pieces = []
    for a, b in base:
        cur = a
        for c, d in holes:
            if d <= cur or c >= b:
                continue
            if c > cur:
                pieces.append((cur, c))
            cur = max(cur, d)
            if cur >= b:
                break
        if cur < b:
            pieces.append((cur, b))
    return pieces


def _overlap(a_start, a_end, b_start, b_end) -> bool:
    """True for positive-duration overlap of two half-open intervals."""
    return max(a_start, b_start) < min(a_end, b_end)


# ----------------------------------------------------------------- windows
def _proven_windows(results: dict[str, PassResult]) -> list[ProvenWindow]:
    windows = []
    for sid, res in results.items():
        for w in res.windows:
            windows.append(ProvenWindow(
                station=sid, start=w.start, end=w.end,
                truncated=w.truncated,
            ))
    return windows


def _components(windows: list[ProvenWindow]):
    """Maximal connected groups under positive-duration overlap.

    Two windows are linked when they overlap over a positive interval and
    belong to different stations (a handoff is only ever between
    stations). Connectivity is computed with a union-find over the
    O(W^2) window pairs; window counts are small (a handful per station
    per day). A component may chain through windows of the same station
    linked by other stations -- every window in it is reachable by a
    valid relay chain, which is exactly what connectivity expresses.
    """
    windows = sorted(windows, key=lambda w: (w.start, w.end, w.station))
    parent = list(range(len(windows)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    for i, wi in enumerate(windows):
        for j in range(i):
            wj = windows[j]
            if wj.station != wi.station and _overlap(
                    wi.start, wi.end, wj.start, wj.end):
                union(i, j)

    groups: dict[int, list[ProvenWindow]] = {}
    for i, w in enumerate(windows):
        groups.setdefault(find(i), []).append(w)
    return [sorted(grp, key=lambda w: (w.start, w.end, w.station))
            for grp in groups.values()]


def _component_span(group: list[ProvenWindow]):
    """Relayable continuous span of a component: [min start, max end).

    On a connected positive-overlap graph every intermediate instant is
    covered by at least one window, and the greedy farthest-reach argument
    guarantees a feasible chain of positive-overlap handoffs from a first
    window to a last one spanning the whole interval.
    """
    return min(w.start for w in group), max(w.end for w in group)


# ------------------------------------------------------------- best chain
def _best_chain(group: list[ProvenWindow], span_start, span_end):
    """Minimum-switch station sequence covering [span_start, span_end).

    Windows are nodes; an edge w_j -> w_i exists when the stations differ
    (every switch is a handoff), w_i positively contains the end of w_j
    (``w_i.start < w_j.end < w_i.end``) and w_i extends the relay strictly
    further in time (``w_i.end > w_j.end``). The incoming window is free
    to have started earlier than the outgoing one. Processing nodes by
    increasing end makes this a DAG shortest path in switches.

    Ties are decided by the lexicographically smallest tuple of station
    IDs of the responsible sequence -- the stable station-ID rule applied
    at every switch, including the choice of the first station.
    Returns (switches, station-id tuple, node sequence) or None.
    """
    nodes = sorted(group, key=lambda w: (w.end, w.start, w.station))

    cost: dict[int, int] = {}
    chain: dict[int, tuple[str, ...]] = {}
    prev: dict[int, Optional[int]] = {}

    for i, wi in enumerate(nodes):
        if wi.start == span_start:
            cost[i] = 0
            chain[i] = (wi.station,)
            prev[i] = None
        for j, wj in enumerate(nodes):
            if j == i or wj.station == wi.station:
                continue
            if cost.get(j) is None:
                continue
            # w_i positively contains w_j's end and reaches further.
            if (wi.start < wj.end < wi.end) and wi.end > wj.end:
                c_new = cost[j] + 1
                ch_new = chain[j] + (wi.station,)
                if cost.get(i) is None or (c_new, ch_new) < (cost[i], chain[i]):
                    cost[i] = c_new
                    chain[i] = ch_new
                    prev[i] = j

    best = None
    for i, wi in enumerate(nodes):
        if wi.end != span_end or cost.get(i) is None:
            continue
        key = (cost[i], chain[i])
        if best is None or key < best[0]:
            best = (key, i)
    if best is None:
        return None

    seq_idx = []
    cur = best[1]
    while cur is not None:
        seq_idx.append(cur)
        cur = prev[cur]
    seq_idx.reverse()
    return cost[best[1]], chain[best[1]], [nodes[i] for i in seq_idx]


# ------------------------------------------------------------- assemble
def _make_plan(group: list[ProvenWindow]) -> Optional[RelayPlan]:
    span_start, span_end = _component_span(group)
    picked = _best_chain(group, span_start, span_end)
    if picked is None:
        return None
    switches, station_tuple, seq = picked

    assignments: list[Assignment] = []
    handoffs: list[Handoff] = []
    cur_start = span_start
    for k, w in enumerate(seq):
        if k + 1 < len(seq):
            # Switch at the instant the outgoing window ends: that instant
            # is still strictly inside the incoming window (positive
            # overlap), so the switch is legal and the tiling gap-free.
            seg_end = w.end
            nxt = seq[k + 1]
            handoffs.append(Handoff(
                time=w.end, from_station=w.station, to_station=nxt.station))
        else:
            seg_end = span_end
        assignments.append(Assignment(
            station=w.station,
            start=cur_start,
            end=seg_end,
            duration_s=(seg_end - cur_start).total_seconds(),
            truncated=w.truncated,
        ))
        cur_start = seg_end

    return RelayPlan(
        start=span_start,
        end=span_end,
        duration_s=(span_end - span_start).total_seconds(),
        switches=switches,
        assignments=assignments,
        handoffs=handoffs,
        stations=station_tuple,
    )


def _intersect(a_list, b_list):
    """Intersection of two interval lists (positive-length pieces only)."""
    out = []
    i = j = 0
    while i < len(a_list) and j < len(b_list):
        lo = max(a_list[i][0], b_list[j][0])
        hi = min(a_list[i][1], b_list[j][1])
        if lo < hi:
            out.append((lo, hi))
        if a_list[i][1] < b_list[j][1]:
            i += 1
        else:
            j += 1
    return out


def _unknown_intervals(search_start, search_end,
                       windows: list[ProvenWindow],
                       results: dict[str, PassResult]) -> list[UnknownInterval]:
    """An instant is *unknown* when no station proves it visible AND at
    least one station registers it in its uncovered_intervals (it could not
    resolve it). That is the intersection of the not-covered remainder with
    the union of every station's unresolved set. Remainder that every
    station proved below threshold is definite non-visibility and is NOT
    listed."""
    covered = _merge([(w.start, w.end) for w in windows])
    gaps = _subtract([(search_start, search_end)], covered)
    if not gaps:
        return []

    sources = []  # (start, end, "station: reason")
    unresolved = []
    for sid, res in results.items():
        for u in res.uncovered:
            a = max(search_start, u.start)
            b = min(search_end, u.end)
            if a < b:
                unresolved.append((a, b))
                sources.append((a, b, f"{sid}: {u.reason}"))
    unresolved = _merge(unresolved)

    unknown = []
    for ga, gb in _intersect(gaps, unresolved):
        hits = sorted({why for a, b, why in sources if a < gb and b > ga})
        unknown.append(UnknownInterval(ga, gb, hits))
    return unknown


def _select_plan(plans: list[RelayPlan]) -> RelayPlan:
    """Longest duration; ties: fewer switches, earlier start, then the
    lexicographically smallest tuple of responsible station IDs."""
    return min(plans, key=lambda p: (-p.duration_s, p.switches,
                                    p.start, p.stations))


# ------------------------------------------------------------------ public
def build_relay(search_start: datetime, search_end: datetime,
                results: dict[str, PassResult]) -> RelayResult:
    """Build the relay plan from per-station :func:`find_passes` results.

    ``results`` maps station ID to its :class:`PassResult` in request
    order. Times must already be normalized to UTC.
    """
    windows = _proven_windows(results)
    unknown = _unknown_intervals(search_start, search_end, windows, results)

    summaries = [
        StationSummary(
            id=sid,
            status=res.status,
            windows=[ProvenWindow(sid, w.start, w.end, w.truncated)
                     for w in res.windows],
            uncovered=list(res.uncovered),
            truncated=any(w.truncated for w in res.windows),
        )
        for sid, res in results.items()
    ]

    plans: list[RelayPlan] = []
    for group in _components(windows):
        plan = _make_plan(group)
        if plan is not None:
            plans.append(plan)

    if plans:
        plan = _select_plan(plans)
        covers_all = (
            plan.start == search_start and plan.end == search_end
            and not unknown
        )
        status = "complete" if covers_all else "partial"
    else:
        plan = None
        status = "failed" if unknown else "complete"

    return RelayResult(
        status=status,
        plan=plan,
        stations=summaries,
        unknown_intervals=unknown,
    )
