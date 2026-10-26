"""Compare already aligned passive-tracer source increments, never raw model files."""

import json
import math
from datetime import datetime, timedelta
from pathlib import Path

from .inventory import sha256

CONTRACT_FIELDS = {
    "tracer_id", "source_id", "emissions_id", "meteorology_case_id",
    "spatial_support_id", "vertical_support_id", "processes_id", "units", "quantity",
}


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a nonempty string.")
    return value


def _utc(value: object) -> datetime:
    text = _text(value, "timestamp")
    parsed = datetime.fromisoformat(text)
    if parsed.utcoffset() != timedelta(0):
        raise ValueError("Timestamps must include a UTC offset (Z or +00:00).")
    return parsed


def load_dataset(path: Path) -> tuple[str, dict[str, str], dict[tuple[datetime, datetime, str], float]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict) or set(data) != {"schema_version", "model", "provenance", "contract", "records"}:
        raise ValueError("Dataset requires schema_version, model, provenance, contract, and records only.")
    if type(data["schema_version"]) is not int or data["schema_version"] != 1:
        raise ValueError("Unsupported dataset schema_version.")
    model = data["model"]
    if model not in ("CMAQ", "AERMOD"):
        raise ValueError("model must be CMAQ or AERMOD.")
    _text(data["provenance"], "provenance")
    raw_contract = data["contract"]
    if not isinstance(raw_contract, dict) or set(raw_contract) != CONTRACT_FIELDS:
        raise ValueError(f"contract must contain exactly {sorted(CONTRACT_FIELDS)}.")
    contract = {key: _text(value, key) for key, value in raw_contract.items()}
    if contract["units"] != "ug/m3" or contract["quantity"] != "passive_tracer_source_increment":
        raise ValueError("Initial comparisons require ug/m3 passive_tracer_source_increment.")
    rows = data["records"]
    if not isinstance(rows, list) or not rows:
        raise ValueError("records must be a nonempty list.")
    values: dict[tuple[datetime, datetime, str], float] = {}
    intervals: dict[str, list[tuple[datetime, datetime]]] = {}
    for row in rows:
        if not isinstance(row, dict) or set(row) != {"start_utc", "end_utc", "support_id", "value"}:
            raise ValueError("Each record requires start_utc, end_utc, support_id, and value only.")
        start, end = _utc(row["start_utc"]), _utc(row["end_utc"])
        support = _text(row["support_id"], "support_id")
        if start >= end:
            raise ValueError("Each averaging interval must have positive duration.")
        value = row["value"]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError("Concentration must be a finite numeric value.")
        key = (start, end, support)
        if key in values:
            raise ValueError(f"Duplicate interval/support record: {key}")
        # Signed source increments are retained; no clipping of subtraction noise.
        values[key] = float(value)
        intervals.setdefault(support, []).append((start, end))
    for support, windows in intervals.items():
        windows.sort()
        if any(second[0] < first[1] for first, second in zip(windows, windows[1:])):
            raise ValueError(f"Overlapping averaging intervals for support {support}.")
    return model, contract, values


def compare(cmaq_path: Path, aermod_path: Path) -> dict[str, object]:
    cmaq_hash, aermod_hash = sha256(cmaq_path), sha256(aermod_path)
    model1, contract1, cmaq = load_dataset(cmaq_path)
    model2, contract2, aermod = load_dataset(aermod_path)
    if (model1, model2) != ("CMAQ", "AERMOD"):
        raise ValueError("Supply CMAQ first and AERMOD second.")
    if contract1 != contract2:
        mismatches = sorted(key for key in CONTRACT_FIELDS if contract1[key] != contract2[key])
        raise ValueError(f"Comparison contracts differ: {', '.join(mismatches)}")
    if cmaq.keys() != aermod.keys():
        raise ValueError("Averaging intervals/support IDs differ; no implicit intersection or interpolation.")
    errors = [cmaq[key] - aermod[key] for key in sorted(cmaq)]
    if not all(math.isfinite(value) for value in errors):
        raise ValueError("Differences overflowed; check concentration scale.")
    count = len(errors)
    # Scale before summing/squaring so large finite inputs do not overflow a sum of squares.
    scale = max(abs(value) for value in errors)
    bias = math.fsum(value / count for value in errors)
    mae = math.fsum(abs(value) / count for value in errors)
    rmse = scale * math.sqrt(math.fsum((value / scale) ** 2 / count for value in errors)) if scale else 0.0
    if cmaq_hash != sha256(cmaq_path) or aermod_hash != sha256(aermod_path):
        raise ValueError("Comparison input changed while being read.")
    return {
        "schema_version": 1, "kind": "aligned_model_comparison", "pairs": count,
        "difference_direction": "CMAQ minus AERMOD", "units": "ug/m3",
        "mean_bias": bias, "mean_absolute_error": mae, "root_mean_square_error": rmse,
        "maximum_absolute_difference": scale, "contract": contract1,
        "inputs": [
            {"model": "CMAQ", "path": str(cmaq_path.resolve()), "sha256": cmaq_hash},
            {"model": "AERMOD", "path": str(aermod_path.resolve()), "sha256": aermod_hash},
        ],
        "scope": "User-declared aligned samples; no physical alignment verification, equivalence, or accuracy ranking.",
    }
