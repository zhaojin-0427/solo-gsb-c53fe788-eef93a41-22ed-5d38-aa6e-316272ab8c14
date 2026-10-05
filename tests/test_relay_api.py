"""Tests for POST /relay-plan and the relay planner."""
from __future__ import annotations

from datetime import datetime, timedelta, timezone

import numpy as np
import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.mechanics import TwoBodyOrbit
from app.geometry import Station
from app.passes import find_passes, PassResult, PassWindow
from app.relay import build_relay
from conftest import EPOCH, circular_leo_state

client = TestClient(app)

T0 = datetime(2025, 1, 1, tzinfo=timezone.utc)


def iso(dt):
    return dt.isoformat().replace("+00:00", "Z")


def _dt(s):
    return T0 + timedelta(seconds=s)


def _win(a, b):
    return PassWindow(
        start=_dt(a), end=_dt(b), peak_time=_dt(a), peak_elevation_deg=20.0,
        peak_azimuth_deg=0.0, rise_azimuth_deg=0.0, set_azimuth_deg=0.0,
        max_range_km=1.0, duration_s=b - a,
        start_open=False, end_open=False, truncated=False,
    )


def _result(windows, status="complete", uncovered=()):
    r = PassResult(status=status)
    r.windows = [_win(a, b) for a, b in windows]
    r.uncovered = list(uncovered)
    return r


def _plan(results, start=0, end=100):
    return build_relay(_dt(start), _dt(end), results)


# --------------------------------------------------------------- HTTP happy
def test_relay_two_stations_single_handoff():
    """Two real overlapping passes: the plan covers the union with exactly
    one handoff inside the positive overlap."""
    epoch, pos, vel = circular_leo_state()
    body = {
        "epoch": iso(epoch),
        "start": "2025-06-01T00:01:40.294555Z",
        "end": "2025-06-01T00:09:41.322998Z",
        "satellite": {"position_km": pos, "velocity_km_s": vel},
        "stations": [
            {"id": "BEIJING", "latitude_deg": 39.9, "longitude_deg": 116.4,
             "altitude_km": 0.05, "min_elevation_deg": 10.0},
            {"id": "TOKYO-BAY", "latitude_deg": 35.5, "longitude_deg": 139.7,
             "altitude_km": 0.02, "min_elevation_deg": 10.0},
        ],
    }
    r = client.post("/relay-plan", json=body)
    assert r.status_code == 200, r.text
    d = r.json()
    assert d["status"] == "complete"
    p = d["plan"]
    assert p["switches"] == 1
    assert [a["station_id"] for a in p["assignments"]] == ["BEIJING", "TOKYO-BAY"]
    # assignments tile [start, end) gap-free
    assert p["assignments"][0]["start"] == d["search_start"]
    assert p["assignments"][-1]["end"] == d["search_end"]
    for a1, a2 in zip(p["assignments"], p["assignments"][1:]):
        assert a1["end"] == a2["start"]
        assert (datetime.fromisoformat(a1["end"].replace("Z", "+00:00"))
                - datetime.fromisoformat(a1["start"].replace("Z", "+00:00"))
                ).total_seconds() > 0
    # the single handoff is at the outgoing window end and strictly inside
    # the incoming window (positive overlap)
    h = p["handoffs"][0]
    assert h["from_station_id"] == "BEIJING"
    assert h["to_station_id"] == "TOKYO-BAY"
    wmap = {s["id"]: s["windows"] for s in d["stations"]}
    assert h["time"] == wmap["BEIJING"][0]["end"]
    assert wmap["TOKYO-BAY"][0]["start"] < h["time"] < wmap["TOKYO-BAY"][0]["end"]
    assert d["unknown_intervals"] == []


def test_relay_wider_window_is_partial_not_complete():
    epoch, pos, vel = circular_leo_state()
    body = {
        "epoch": iso(epoch),
        "start": iso(EPOCH),
        "end": iso(EPOCH + timedelta(minutes=30)),
        "satellite": {"position_km": pos, "velocity_km_s": vel},
        "stations": [
            {"id": "BEIJING", "latitude_deg": 39.9, "longitude_deg": 116.4,
             "altitude_km": 0.05, "min_elevation_deg": 10.0},
            {"id": "TOKYO-BAY", "latitude_deg": 35.5, "longitude_deg": 139.7,
             "altitude_km": 0.02, "min_elevation_deg": 10.0},
        ],
    }
    d = client.post("/relay-plan", json=body).json()
    # definite invisible gaps exist and are NOT listed as unknown
    assert d["status"] == "partial"
    assert d["unknown_intervals"] == []
    assert d["plan"] is not None
    assert d["plan"]["duration_s"] < 1800.0


def test_relay_no_visibility_empty_plan_complete():
    """Everybody proven below threshold: empty plan, status complete, no
    unknown intervals (must not claim failed or full coverage)."""
    epoch, pos, vel = circular_leo_state()
    body = {
        "epoch": iso(epoch),
        "start": iso(EPOCH),
        "end": iso(EPOCH + timedelta(hours=1)),
        "satellite": {"position_km": pos, "velocity_km_s": vel},
        "stations": [
            {"id": "S1", "latitude_deg": 39.9, "longitude_deg": 116.4,
             "altitude_km": 0.05, "min_elevation_deg": 89.0},
            {"id": "S2", "latitude_deg": 35.5, "longitude_deg": 139.7,
             "altitude_km": 0.02, "min_elevation_deg": 89.0},
        ],
    }
    d = client.post("/relay-plan", json=body).json()
    assert d["status"] == "complete"
    assert d["plan"] is None
    assert d["unknown_intervals"] == []
    assert all(s["windows"] == [] for s in d["stations"])


def test_relay_grazing_uncovered_is_unknown_not_invisible():
    """One station grazes (unresolved), the other is proven nowhere: the
    gap is unknown and the plan is empty/failed -- never complete."""
    epoch, pos, vel = circular_leo_state()
    orb = TwoBodyOrbit(np.array(pos), np.array(vel))
    res = find_passes(orb, Station(39.9, 116.4, 0.05), epoch,
                      EPOCH, EPOCH + timedelta(hours=1), 10.0)
    w0 = res.windows[0]
    peak, tp = w0.peak_elevation_deg, w0.peak_time
    body = {
        "epoch": iso(epoch),
        "start": iso(tp - timedelta(seconds=30)),
        "end": iso(tp + timedelta(seconds=30)),
        "satellite": {"position_km": pos, "velocity_km_s": vel},
        "stations": [
            {"id": "GRAZE", "latitude_deg": 39.9, "longitude_deg": 116.4,
             "altitude_km": 0.05, "min_elevation_deg": peak + 1e-7},
            {"id": "FAR", "latitude_deg": 0.0, "longitude_deg": 0.0,
             "altitude_km": 0.0, "min_elevation_deg": 80.0},
        ],
    }
    d = client.post("/relay-plan", json=body).json()
    assert d["plan"] is None
    assert d["status"] == "failed"
    assert len(d["unknown_intervals"]) == 1
    u = d["unknown_intervals"][0]
    assert u["start"] < u["end"]
    assert any("GRAZE" in reason for reason in u["reasons"])


def test_relay_one_uncovered_other_visible_still_relays():
    """The grazing station cannot resolve an interval, but the same interval
    is proven visible by a normal station: it is covered (not unknown) and
    the plan proceeds on the visible station."""
    epoch, pos, vel = circular_leo_state()
    orb = TwoBodyOrbit(np.array(pos), np.array(vel))
    res = find_passes(orb, Station(39.9, 116.4, 0.05), epoch,
                      EPOCH, EPOCH + timedelta(hours=1), 10.0)
    w0 = res.windows[0]
    peak, tp = w0.peak_elevation_deg, w0.peak_time
    body = {
        "epoch": iso(epoch),
        "start": iso(tp - timedelta(seconds=30)),
        "end": iso(tp + timedelta(seconds=30)),
        "satellite": {"position_km": pos, "velocity_km_s": vel},
        "stations": [
            {"id": "GRAZE", "latitude_deg": 39.9, "longitude_deg": 116.4,
             "altitude_km": 0.05, "min_elevation_deg": peak + 1e-7},
            {"id": "NORMAL", "latitude_deg": 39.9, "longitude_deg": 116.4,
             "altitude_km": 0.05, "min_elevation_deg": 10.0},
        ],
    }
    d = client.post("/relay-plan", json=body).json()
    assert d["status"] == "complete"
    assert d["unknown_intervals"] == []
    assert d["plan"]["stations"] == ["NORMAL"]
    assert d["plan"]["duration_s"] > 0


# --------------------------------------------------------------- HTTP 422
def test_relay_duplicate_ids_422():
    epoch, pos, vel = circular_leo_state()
    body = {
        "epoch": iso(epoch), "start": iso(EPOCH),
        "end": iso(EPOCH + timedelta(hours=1)),
        "satellite": {"position_km": pos, "velocity_km_s": vel},
        "stations": [
            {"id": "X", "latitude_deg": 39.9, "longitude_deg": 116.4,
             "altitude_km": 0.05, "min_elevation_deg": 10.0},
            {"id": "X", "latitude_deg": 35.5, "longitude_deg": 139.7,
             "altitude_km": 0.02, "min_elevation_deg": 10.0},
        ],
    }
    r = client.post("/relay-plan", json=body)
    assert r.status_code == 422
    assert "unique" in r.json()["detail"][0]["msg"]


def test_relay_empty_stations_422():
    epoch, pos, vel = circular_leo_state()
    body = {
        "epoch": iso(epoch), "start": iso(EPOCH),
        "end": iso(EPOCH + timedelta(hours=1)),
        "satellite": {"position_km": pos, "velocity_km_s": vel},
        "stations": [],
    }
    assert client.post("/relay-plan", json=body).status_code == 422


@pytest.mark.parametrize("bad_id", ["", "   ", "\t"])
def test_relay_blank_id_422(bad_id):
    epoch, pos, vel = circular_leo_state()
    body = {
        "epoch": iso(epoch), "start": iso(EPOCH),
        "end": iso(EPOCH + timedelta(hours=1)),
        "satellite": {"position_km": pos, "velocity_km_s": vel},
        "stations": [
            {"id": bad_id, "latitude_deg": 39.9, "longitude_deg": 116.4,
             "altitude_km": 0.05, "min_elevation_deg": 10.0}],
    }
    assert client.post("/relay-plan", json=body).status_code == 422


def test_relay_degenerate_orbit_422():
    epoch, pos, vel = circular_leo_state()
    body = {
        "epoch": iso(epoch), "start": iso(EPOCH),
        "end": iso(EPOCH + timedelta(hours=1)),
        "satellite": {"position_km": pos, "velocity_km_s": [x * 1000 for x in vel]},
        "stations": [
            {"id": "X", "latitude_deg": 39.9, "longitude_deg": 116.4,
             "altitude_km": 0.05, "min_elevation_deg": 10.0}],
    }
    r = client.post("/relay-plan", json=body)
    assert r.status_code == 422
    assert r.json()["detail"][0]["type"] == "orbit_degenerate_or_invalid"


def test_relay_impact_orbit_422():
    # Circular-speed state at a radius 50 km below the surface: a bound
    # orbit whose (circular) perigee lies inside the Earth.
    from app.mechanics import MU_EARTH_KM, EARTH_RADIUS_KM
    r = EARTH_RADIUS_KM - 50.0
    vcirc = (MU_EARTH_KM / r) ** 0.5
    body = {
        "epoch": iso(EPOCH), "start": iso(EPOCH),
        "end": iso(EPOCH + timedelta(hours=1)),
        "satellite": {"position_km": [r, 0.0, 0.0],
                      "velocity_km_s": [0.0, vcirc, 0.0]},
        "stations": [
            {"id": "X", "latitude_deg": 39.9, "longitude_deg": 116.4,
             "altitude_km": 0.05, "min_elevation_deg": 10.0}],
    }
    r2 = client.post("/relay-plan", json=body)
    assert r2.status_code == 422
    assert r2.json()["detail"][0]["type"] == "orbit_impact_earth"


def test_relay_span_over_seven_days_422():
    epoch, pos, vel = circular_leo_state()
    body = {
        "epoch": iso(epoch), "start": iso(EPOCH),
        "end": iso(EPOCH + timedelta(days=7, seconds=1)),
        "satellite": {"position_km": pos, "velocity_km_s": vel},
        "stations": [
            {"id": "X", "latitude_deg": 39.9, "longitude_deg": 116.4,
             "altitude_km": 0.05, "min_elevation_deg": 10.0}],
    }
    assert client.post("/relay-plan", json=body).status_code == 422


def test_relay_naive_datetime_422():
    epoch, pos, vel = circular_leo_state()
    body = {
        "epoch": "2025-06-01T00:00:00", "start": iso(EPOCH),
        "end": iso(EPOCH + timedelta(hours=1)),
        "satellite": {"position_km": pos, "velocity_km_s": vel},
        "stations": [
            {"id": "X", "latitude_deg": 39.9, "longitude_deg": 116.4,
             "altitude_km": 0.05, "min_elevation_deg": 10.0}],
    }
    assert client.post("/relay-plan", json=body).status_code == 422


@pytest.mark.parametrize("field,value", [
    ("latitude_deg", 90.1), ("longitude_deg", -181.0),
    ("altitude_km", -13.0), ("min_elevation_deg", 90.0),
])
def test_relay_bad_station_range_422(field, value):
    epoch, pos, vel = circular_leo_state()
    body = {
        "epoch": iso(epoch), "start": iso(EPOCH),
        "end": iso(EPOCH + timedelta(hours=1)),
        "satellite": {"position_km": pos, "velocity_km_s": vel},
        "stations": [
            {"id": "X", "latitude_deg": 39.9, "longitude_deg": 116.4,
             "altitude_km": 0.05, "min_elevation_deg": 10.0}],
    }
    body["stations"][0][field] = value
    assert client.post("/relay-plan", json=body).status_code == 422


def test_relay_extra_field_422():
    epoch, pos, vel = circular_leo_state()
    body = {
        "epoch": iso(epoch), "start": iso(EPOCH),
        "end": iso(EPOCH + timedelta(hours=1)),
        "satellite": {"position_km": pos, "velocity_km_s": vel},
        "stations": [
            {"id": "X", "latitude_deg": 39.9, "longitude_deg": 116.4,
             "altitude_km": 0.05, "min_elevation_deg": 10.0}],
        "bogus": 1,
    }
    assert client.post("/relay-plan", json=body).status_code == 422


# ------------------------------------------------------- planner unit rules
def test_touching_windows_cannot_hand_off():
    """Half-open windows touching at one instant give no positive overlap,
    so they cannot be one continuous relay."""
    rr = _plan({"a": _result([(0, 3)]), "b": _result([(3, 6)])}, 0, 6)
    assert rr.plan.duration_s == 3.0
    assert rr.plan.switches == 0


def test_positive_overlap_handoff():
    rr = _plan({"a": _result([(0, 4)]), "b": _result([(2, 6)])}, 0, 6)
    assert rr.plan.duration_s == 6.0
    assert rr.plan.switches == 1
    assert rr.plan.stations == ("a", "b")
    h = rr.plan.handoffs[0]
    assert h.time == _dt(4)
    assert h.from_station == "a" and h.to_station == "b"


def test_switch_back_through_same_station():
    """a1 -> b -> a2 covers a span longer than any single window."""
    rr = _plan({"x": _result([(0, 5), (4, 10)]), "y": _result([(4.5, 9)])},
               0, 10)
    assert rr.plan.stations == ("x", "y", "x")
    assert rr.plan.switches == 2
    assert rr.plan.duration_s == 10.0


def test_tie_fewer_switches_wins():
    # 1-switch chain a->b reaches 6; c offers only a 2-switch alternative.
    rr = _plan({"a": _result([(0, 3)]), "b": _result([(2, 6)]),
                "c": _result([(2, 5)])}, 0, 6)
    assert rr.plan.switches == 1
    assert rr.plan.stations == ("a", "b")


def test_tie_earlier_start_then_station_id():
    # two disjoint relayable components of equal length and switches
    rr = _plan({"a": _result([(0, 4)]), "b": _result([(2, 6)]),
                "c": _result([(10, 14)]), "d": _result([(12, 16)])}, 0, 20)
    assert rr.plan.start == _dt(0)

    # same span & switches: lexicographically smallest first station wins
    rr2 = _plan({"z": _result([(0, 4)]), "a": _result([(0, 4)]),
                 "b": _result([(3, 6)])}, 0, 6)
    assert rr2.plan.stations[0] == "a"


def test_longest_component_is_chosen():
    rr = _plan({"a": _result([(0, 4)]), "b": _result([(2, 6)]),
                "c": _result([(10, 18)]), "d": _result([(11, 20)])}, 0, 20)
    assert rr.plan.duration_s == 10.0
    assert rr.plan.start == _dt(10)


def test_unknown_interval_listed_only_where_unresolved():
    from app.passes import UncoveredInterval
    unc = [UncoveredInterval(_dt(7), _dt(9), "grazing")]
    # a visible window covers 0..6; 6..7 proven invisible; 7..9 unresolved
    rr = _plan({"a": _result([(0, 6)], uncovered=unc)}, 0, 10)
    assert rr.status == "partial"
    assert rr.plan.start == _dt(0) and rr.plan.end == _dt(6)
    assert len(rr.unknown_intervals) == 1
    u = rr.unknown_intervals[0]
    assert u.start == _dt(7) and u.end == _dt(9)
    # the 6..7 definite-invisible gap is not unknown
    assert all(not (u.start <= _dt(6) < u.end) for u in rr.unknown_intervals)


def test_uncovered_covered_by_other_station_is_not_unknown():
    from app.passes import UncoveredInterval
    unc = [UncoveredInterval(_dt(2), _dt(4), "grazing")]
    rr = _plan({
        "g": _result([], uncovered=unc),
        "ok": _result([(0, 6)]),
    }, 0, 6)
    assert rr.status == "complete"
    assert rr.unknown_intervals == []
    assert rr.plan.stations == ("ok",)


def test_all_invisible_with_uncovered_is_failed():
    from app.passes import UncoveredInterval
    unc = [UncoveredInterval(_dt(0), _dt(2), "grazing")]
    rr = _plan({"g": _result([], uncovered=unc)}, 0, 5)
    assert rr.plan is None
    assert rr.status == "failed"
    assert len(rr.unknown_intervals) == 1
