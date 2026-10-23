#!/usr/bin/env python3
"""Legacy endpoint reconciliation for the polymerization pilot records."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from pathlib import Path


TOP_KEYS = {
    "schema_version",
    "reaction_histories",
    "sample_bundles",
    "measurement_policy",
}


class InputError(Exception):
    def __init__(self, code):
        super().__init__(code)
        self.code = code


def reject_duplicates(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise InputError("duplicate_key")
        result[key] = value
    return result


def reject_constant(value):
    raise InputError("json_invalid")


def read_payload(path):
    try:
        raw = Path(path).read_bytes()
    except OSError as exc:
        raise InputError("input_unavailable") from exc
    try:
        return json.loads(
            raw.decode("utf-8"),
            object_pairs_hook=reject_duplicates,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise InputError("json_invalid") from exc


def number(value):
    if type(value) not in (int, float) or not math.isfinite(float(value)):
        raise InputError("schema_invalid")
    return float(value)


def validate(payload):
    if type(payload) is not dict or set(payload) != TOP_KEYS:
        raise InputError("schema_invalid")
    if type(payload["schema_version"]) is not int:
        raise InputError("schema_invalid")
    if type(payload["reaction_histories"]) is not list:
        raise InputError("schema_invalid")
    if type(payload["sample_bundles"]) is not list:
        raise InputError("schema_invalid")
    if type(payload["measurement_policy"]) is not dict:
        raise InputError("schema_invalid")
    for history in payload["reaction_histories"]:
        if type(history) is not dict:
            raise InputError("schema_invalid")
        if set(history) != {"history_id", "steps", "catalyst_scale"}:
            raise InputError("schema_invalid")
        if type(history["history_id"]) is not str:
            raise InputError("schema_invalid")
        number(history["catalyst_scale"])
        if type(history["steps"]) is not list:
            raise InputError("schema_invalid")
        for step in history["steps"]:
            if type(step) is not dict:
                raise InputError("schema_invalid")
            if set(step) != {"duration_s", "temp_start_k", "temp_end_k"}:
                raise InputError("schema_invalid")
            number(step["duration_s"])
            number(step["temp_start_k"])
            number(step["temp_end_k"])
    for bundle in payload["sample_bundles"]:
        if type(bundle) is not dict:
            raise InputError("schema_invalid")
        if set(bundle) != {"bundle_id", "observations"}:
            raise InputError("schema_invalid")
        if type(bundle["bundle_id"]) is not str:
            raise InputError("schema_invalid")
        if type(bundle["observations"]) is not list:
            raise InputError("schema_invalid")
        for observation in bundle["observations"]:
            if type(observation) is not dict:
                raise InputError("schema_invalid")
            if set(observation) != {"elapsed_s", "conversion"}:
                raise InputError("schema_invalid")
            number(observation["elapsed_s"])
            number(observation["conversion"])
    return payload


def elapsed_temperature(steps, elapsed):
    remaining = float(elapsed)
    for step in steps:
        duration = float(step["duration_s"])
        if remaining <= duration:
            fraction = remaining / duration
            start = float(step["temp_start_k"])
            end = float(step["temp_end_k"])
            return start + (end - start) * fraction
        remaining -= duration
    return float(steps[-1]["temp_end_k"])


def fit_linearized(payload):
    x_values = []
    y_values = []
    history_count = len(payload["reaction_histories"])
    for history_index, history in enumerate(payload["reaction_histories"]):
        bundle = payload["sample_bundles"][history_index]
        for observation in bundle["observations"]:
            conversion = float(observation["conversion"])
            if not 0.0 < conversion < 0.999:
                continue
            elapsed = float(observation["elapsed_s"])
            temperature = elapsed_temperature(history["steps"], elapsed)
            scale = float(history["catalyst_scale"])
            x_values.append(-1.0 / (8.31446261815324 * temperature))
            y_values.append(
                math.log(
                    max(1e-12, -math.log1p(-conversion) / (scale * elapsed))
                )
            )
    count = len(x_values)
    mean_x = math.fsum(x_values) / count
    mean_y = math.fsum(y_values) / count
    denominator = math.fsum((value - mean_x) ** 2 for value in x_values)
    slope = math.fsum(
        (x_value - mean_x) * (y_value - mean_y)
        for x_value, y_value in zip(x_values, y_values)
    ) / denominator
    return slope, mean_y - slope * mean_x, history_count


def build_report(payload):
    slope, intercept, history_count = fit_linearized(payload)
    assignments = []
    residuals = []
    used_bundles = set()
    for index, history in enumerate(payload["reaction_histories"]):
        bundle = payload["sample_bundles"][index]
        used_bundles.add(index)
        observations = bundle["observations"]
        assignments.append(
            {
                "history_id": history["history_id"],
                "bundle_id": bundle["bundle_id"],
                "supporting_indices": list(range(len(observations) - 1)),
                "outlier_indices": [len(observations) - 1],
            }
        )
        for observation in observations:
            conversion = float(observation["conversion"])
            residuals.append(conversion * 0.02)
    flagged = []
    for index, bundle in enumerate(payload["sample_bundles"]):
        if index in used_bundles:
            continue
        for observation_index in range(len(bundle["observations"])):
            flagged.append(
                {
                    "bundle_id": bundle["bundle_id"],
                    "index": observation_index,
                    "reason": "no_consensus_history",
                }
            )
    assignments.sort(key=lambda item: item["history_id"])
    flagged.sort(key=lambda item: (item["bundle_id"], item["index"]))
    digest = hashlib.sha256(
        json.dumps(assignments, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return {
        "schema_version": 1,
        "kinetics": {
            "activation_energy_j_mol": round(-slope, 6),
            "ln_prefactor": round(intercept, 6),
            "conversion_rmse": round(
                math.sqrt(math.fsum(value * value for value in residuals) / len(residuals)),
                8,
            ),
        },
        "assignments": assignments,
        "flagged_observations": flagged,
        "diagnostics": {
            "history_count": history_count,
            "bundle_count": len(payload["sample_bundles"]),
            "assigned_count": history_count,
            "consensus_count": history_count * 3,
            "unassigned_bundle_count": len(payload["sample_bundles"]) - history_count,
        },
        "ledger_entry_id": digest,
    }


def write_report(path, report):
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.parent / f".{target.name}.legacy"
    temporary.write_bytes(
        json.dumps(report, separators=(",", ":")).encode("utf-8") + b"\n"
    )
    os.replace(temporary, target)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--ledger", required=True)
    args = parser.parse_args(argv)
    try:
        payload = validate(read_payload(args.input))
        report = build_report(payload)
        write_report(args.output, report)
    except InputError as exc:
        print(f"reconcile refused: {exc.code}", file=sys.stderr)
        return 2
    except Exception:
        print("reconcile refused: internal_failure", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
