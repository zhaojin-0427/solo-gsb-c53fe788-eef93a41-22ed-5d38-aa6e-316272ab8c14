"""FastAPI application for satellite pass prediction."""
from __future__ import annotations

import math

import numpy as np
from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse

from .geometry import Station
from .mechanics import EARTH_RADIUS_KM, OrbitError, TwoBodyOrbit
from .passes import find_passes
from .relay import build_relay
from .schemas import (
    PassRequest,
    PassResponse,
    RelayRequest,
    RelayResponse,
)

app = FastAPI(
    title="Satellite Pass Observation API",
    version="1.0.0",
    description=(
        "Predict visibility windows of an Earth satellite from a ground "
        "station. Two-body Kepler propagation, fixed WGS84 Earth rotation, "
        "adaptive stepping with bisection-level crossing accuracy."
    ),
)


@app.exception_handler(OrbitError)
async def orbit_error_handler(request: Request, exc: OrbitError):
    return JSONResponse(
        status_code=422,
        content={
            "detail": [
                {
                    "type": "orbit_degenerate_or_invalid",
                    "loc": ["body", "satellite"],
                    "msg": str(exc),
                }
            ]
        },
    )


def _build_orbit(pos, vel) -> TwoBodyOrbit:
    return TwoBodyOrbit(np.array(pos, dtype=float), np.array(vel, dtype=float))


def _orbit_impact_response(orbit: TwoBodyOrbit) -> JSONResponse:
    return JSONResponse(
        status_code=422,
        content={
            "detail": [
                {
                    "type": "orbit_impact_earth",
                    "loc": ["body", "satellite"],
                    "msg": (
                        f"orbit perigee radius {orbit.periapsis:.3f} km is "
                        f"below the Earth surface (~{EARTH_RADIUS_KM:.3f} km); "
                        "two-body propagation is not physically meaningful"
                    ),
                }
            ]
        },
    )


@app.get("/health")
def health():
    return {"status": "ok"}


@app.post("/passes", response_model=PassResponse)
def passes(req: PassRequest):
    orbit = _build_orbit(req.satellite.position_km,
                         req.satellite.velocity_km_s)

    # Physical reality check: an ellipse whose perigee is inside the Earth
    # cannot be propagated meaningfully for observations.
    min_surface = EARTH_RADIUS_KM - 14.0  # account for deepest ocean trench
    if orbit.periapsis < min_surface:
        return _orbit_impact_response(orbit)

    station = Station(
        lat_deg=req.station.latitude_deg,
        lon_deg=req.station.longitude_deg,
        alt_km=req.station.altitude_km,
    )

    result = find_passes(
        orbit=orbit,
        station=station,
        epoch=req.epoch,
        t_start=req.start,
        t_end=req.end,
        min_elevation_deg=req.min_elevation_deg,
    )

    span_s = (req.end - req.start).total_seconds()
    uncovered_s = sum(
        (u.end - u.start).total_seconds() for u in result.uncovered
    )
    diagnostics = {
        "orbit": {
            "semi_major_axis_km": orbit.a,
            "eccentricity": orbit.e,
            "period_s": orbit.period,
            "perigee_radius_km": orbit.periapsis,
            "apogee_radius_km": orbit.apoapsis,
            "perigee_altitude_km": orbit.periapsis - EARTH_RADIUS_KM,
            "apogee_altitude_km": orbit.apoapsis - EARTH_RADIUS_KM,
        },
        "coarse_samples": result.coarse_samples,
        "dense_samples": result.dense_samples,
        "search_span_s": span_s,
        "uncovered_seconds": uncovered_s,
        "coverage_fraction": max(0.0, 1.0 - uncovered_s / span_s)
        if span_s > 0 else 1.0,
        "notes": (
            "Elevation crossings refined by bisection; peaks by golden-section "
            "search. Half-open search interval [start, end). Below-threshold "
            "intervals count as covered only when a Lipschitz rate bound "
            "proves them invisible; unprovable brackets are bisected down to "
            "the 1 us reporting resolution to expose hidden crossings, and "
            "whatever stays unresolvable (grazing contacts) is reported in "
            "uncovered_intervals instead of being silently skipped."
        ),
    }

    return PassResponse(
        status=result.status,
        count=len(result.windows),
        windows=[
            {
                "start": w.start,
                "end": w.end,
                "peak_time": w.peak_time,
                "peak_elevation_deg": w.peak_elevation_deg,
                "peak_azimuth_deg": _nan_none(w.peak_azimuth_deg),
                "rise_azimuth_deg": _nan_none(w.rise_azimuth_deg),
                "set_azimuth_deg": _nan_none(w.set_azimuth_deg),
                "max_range_km": w.max_range_km,
                "duration_s": w.duration_s,
                "start_open": w.start_open,
                "end_open": w.end_open,
                "truncated": w.truncated,
            }
            for w in result.windows
        ],
        uncovered_intervals=[
            {"start": u.start, "end": u.end, "reason": u.reason}
            for u in result.uncovered
        ],
        search_start=req.start,
        search_end=req.end,
        min_elevation_deg=req.min_elevation_deg,
        diagnostics=diagnostics,
    )


def _nan_none(x: float):
    return None if x is None or (isinstance(x, float) and math.isnan(x)) else x


@app.post("/relay-plan", response_model=RelayResponse)
def relay_plan(req: RelayRequest):
    """Plan the longest continuous proven-visibility relay over the
    stations, reusing the same two-body orbit and find_passes machinery
    as POST /passes (once per station)."""
    orbit = _build_orbit(req.satellite.position_km,
                         req.satellite.velocity_km_s)

    # Same physical reality check as /passes: a sub-surface perigee cannot
    # be propagated meaningfully for any station.
    min_surface = EARTH_RADIUS_KM - 14.0
    if orbit.periapsis < min_surface:
        return _orbit_impact_response(orbit)

    results = {}
    for st in req.stations:
        station = Station(
            lat_deg=st.latitude_deg,
            lon_deg=st.longitude_deg,
            alt_km=st.altitude_km,
        )
        results[st.id] = find_passes(
            orbit=orbit,
            station=station,
            epoch=req.epoch,
            t_start=req.start,
            t_end=req.end,
            min_elevation_deg=st.min_elevation_deg,
        )

    relay = build_relay(req.start, req.end, results)

    plan_out = None
    if relay.plan is not None:
        p = relay.plan
        plan_out = {
            "start": p.start,
            "end": p.end,
            "duration_s": p.duration_s,
            "switches": p.switches,
            "assignments": [
                {
                    "station_id": a.station,
                    "start": a.start,
                    "end": a.end,
                    "duration_s": a.duration_s,
                    "truncated": a.truncated,
                }
                for a in p.assignments
            ],
            "handoffs": [
                {
                    "time": h.time,
                    "from_station_id": h.from_station,
                    "to_station_id": h.to_station,
                }
                for h in p.handoffs
            ],
            "stations": list(p.stations),
        }

    return RelayResponse(
        status=relay.status,
        search_start=req.start,
        search_end=req.end,
        plan=plan_out,
        stations=[
            {
                "id": s.id,
                "status": s.status,
                "truncated": s.truncated,
                "windows": [
                    {"start": w.start, "end": w.end, "truncated": w.truncated}
                    for w in s.windows
                ],
                "uncovered_intervals": [
                    {"start": u.start, "end": u.end, "reason": u.reason}
                    for u in s.uncovered
                ],
            }
            for s in relay.stations
        ],
        unknown_intervals=[
            {"start": u.start, "end": u.end, "reasons": u.reasons}
            for u in relay.unknown_intervals
        ],
    )
