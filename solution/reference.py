#!/usr/bin/env python3
"""Reconcile unlabeled polymerization samples with thermal histories."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import tempfile
from pathlib import Path


GAS_CONSTANT = 8.31446261815324
GAUSS_NODES = (
    -0.9739065285171717,
    -0.8650633666889845,
    -0.6794095682990244,
    -0.4333953941292472,
    -0.1488743389816312,
    0.1488743389816312,
    0.4333953941292472,
    0.6794095682990244,
    0.8650633666889845,
    0.9739065285171717,
)
GAUSS_WEIGHTS = (
    0.0666713443086881,
    0.1494513491505806,
    0.2190863625159820,
    0.2692667193099963,
    0.2955242247147529,
    0.2955242247147529,
    0.2692667193099963,
    0.2190863625159820,
    0.1494513491505806,
    0.0666713443086881,
)
TOP_KEYS = {
    "schema_version",
    "reaction_histories",
    "sample_bundles",
    "measurement_policy",
}
HISTORY_KEYS = {"history_id", "steps", "catalyst_scale"}
STEP_KEYS = {"duration_s", "temp_start_k", "temp_end_k"}
BUNDLE_KEYS = {"bundle_id", "observations"}
OBSERVATION_KEYS = {"elapsed_s", "conversion"}
POLICY_KEYS = {
    "activation_energy_bounds_j_mol",
    "ln_prefactor_bounds",
    "conversion_floor",
    "consensus_min_points",
    "conversion_tolerance",
    "max_assignments",
}
LEDGER_KEYS = {
    "schema_version",
    "entry_id",
    "activation_energy_j_mol",
    "ln_prefactor",
    "assignments_sha256",
}


class Refusal(Exception):
    """A public input or publication failure with a stable code."""

    def __init__(self, code):
        super().__init__(code)
        self.code = code


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise Refusal("duplicate_key")
        result[key] = value
    return result


def _parse_constant(value):
    raise Refusal("json_invalid")


def read_json(path):
    try:
        raw = Path(path).read_bytes()
    except OSError as exc:
        raise Refusal("input_unavailable") from exc
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise Refusal("json_invalid") from exc
    try:
        return json.loads(
            text,
            object_pairs_hook=_pairs,
            parse_constant=_parse_constant,
        )
    except (json.JSONDecodeError, ValueError) as exc:
        if isinstance(exc, Refusal):
            raise
        raise Refusal("json_invalid") from exc


def require_object(value, keys):
    if type(value) is not dict or set(value) != keys:
        raise Refusal("schema_invalid")


def require_number(value):
    if type(value) not in (int, float) or not math.isfinite(float(value)):
        raise Refusal("schema_invalid")
    return float(value)


def require_int(value):
    if type(value) is not int:
        raise Refusal("schema_invalid")
    return value


def require_string(value):
    if type(value) is not str or not value or value != value.strip():
        raise Refusal("schema_invalid")
    return value


def validate_dataset(payload):
    require_object(payload, TOP_KEYS)
    if require_int(payload["schema_version"]) != 1:
        raise Refusal("schema_invalid")

    histories = payload["reaction_histories"]
    bundles = payload["sample_bundles"]
    policy = payload["measurement_policy"]
    if type(histories) is not list or not 8 <= len(histories) <= 30:
        raise Refusal("constraint_invalid")
    if type(bundles) is not list or not 8 <= len(bundles) <= 32:
        raise Refusal("constraint_invalid")
    require_object(policy, POLICY_KEYS)

    history_ids = []
    for history in histories:
        require_object(history, HISTORY_KEYS)
        history_id = require_string(history["history_id"])
        if len(history_id) < 6 or len(history_id) > 24:
            raise Refusal("constraint_invalid")
        history_ids.append(history_id)
        scale = require_number(history["catalyst_scale"])
        if not 0.2 <= scale <= 3.0:
            raise Refusal("constraint_invalid")
        steps = history["steps"]
        if type(steps) is not list or not 2 <= len(steps) <= 8:
            raise Refusal("constraint_invalid")
        for step in steps:
            require_object(step, STEP_KEYS)
            duration = require_number(step["duration_s"])
            low = require_number(step["temp_start_k"])
            high = require_number(step["temp_end_k"])
            if duration <= 0.0 or duration > 7200.0:
                raise Refusal("constraint_invalid")
            if not (280.0 <= low <= 500.0):
                raise Refusal("constraint_invalid")
            if not (280.0 <= high <= 500.0):
                raise Refusal("constraint_invalid")
    if len(set(history_ids)) != len(history_ids):
        raise Refusal("constraint_invalid")

    bundle_ids = []
    for bundle in bundles:
        require_object(bundle, BUNDLE_KEYS)
        bundle_id = require_string(bundle["bundle_id"])
        if len(bundle_id) < 6 or len(bundle_id) > 24:
            raise Refusal("constraint_invalid")
        bundle_ids.append(bundle_id)
        observations = bundle["observations"]
        if type(observations) is not list or not 4 <= len(observations) <= 6:
            raise Refusal("constraint_invalid")
        elapsed_values = []
        for observation in observations:
            require_object(observation, OBSERVATION_KEYS)
            elapsed = require_number(observation["elapsed_s"])
            conversion = require_number(observation["conversion"])
            if elapsed <= 0.0 or not 0.0 <= conversion < 1.0:
                raise Refusal("constraint_invalid")
            elapsed_values.append(elapsed)
        if len(set(elapsed_values)) != len(elapsed_values):
            raise Refusal("constraint_invalid")
    if len(set(bundle_ids)) != len(bundle_ids):
        raise Refusal("constraint_invalid")

    ea_bounds = policy["activation_energy_bounds_j_mol"]
    ln_bounds = policy["ln_prefactor_bounds"]
    if type(ea_bounds) is not list or len(ea_bounds) != 2:
        raise Refusal("schema_invalid")
    if type(ln_bounds) is not list or len(ln_bounds) != 2:
        raise Refusal("schema_invalid")
    ea_low, ea_high = (require_number(value) for value in ea_bounds)
    ln_low, ln_high = (require_number(value) for value in ln_bounds)
    floor = require_number(policy["conversion_floor"])
    tolerance = require_number(policy["conversion_tolerance"])
    min_points = require_int(policy["consensus_min_points"])
    max_assignments = require_int(policy["max_assignments"])
    if not (30000.0 <= ea_low < ea_high <= 180000.0):
        raise Refusal("constraint_invalid")
    if not (10.0 <= ln_low < ln_high <= 35.0):
        raise Refusal("constraint_invalid")
    if not (0.0 < floor < 0.01):
        raise Refusal("constraint_invalid")
    if not (0.000001 <= tolerance <= 0.02):
        raise Refusal("constraint_invalid")
    if min_points < 3:
        raise Refusal("constraint_invalid")
    if any(len(bundle["observations"]) - 1 < min_points for bundle in bundles):
        raise Refusal("constraint_invalid")
    if max_assignments != len(histories):
        raise Refusal("constraint_invalid")
    if max_assignments > 30:
        raise Refusal("constraint_invalid")
    return payload


def thermal_dose(steps, elapsed, activation_energy):
    remaining = float(elapsed)
    total = 0.0
    for step in steps:
        if remaining <= 0.0:
            break
        duration = float(step["duration_s"])
        used = min(remaining, duration)
        fraction = used / duration
        start = float(step["temp_start_k"])
        end = start + (float(step["temp_end_k"]) - start) * fraction
        for node, weight in zip(GAUSS_NODES, GAUSS_WEIGHTS):
            temperature = start + (end - start) * 0.5 * (node + 1.0)
            total += 0.5 * used * weight * math.exp(
                -activation_energy / (GAS_CONSTANT * temperature)
            )
        remaining -= used
    return total


def hungarian(cost):
    """Return a minimum-cost column assignment for a square matrix."""
    n = len(cost)
    u = [0.0] * (n + 1)
    v = [0.0] * (n + 1)
    p = [0] * (n + 1)
    way = [0] * (n + 1)
    for row in range(1, n + 1):
        p[0] = row
        column0 = 0
        minv = [math.inf] * (n + 1)
        used = [False] * (n + 1)
        while True:
            used[column0] = True
            row0 = p[column0]
            delta = math.inf
            column1 = 0
            for column in range(1, n + 1):
                if used[column]:
                    continue
                current = cost[row0 - 1][column - 1] - u[row0] - v[column]
                if current < minv[column]:
                    minv[column] = current
                    way[column] = column0
                if minv[column] < delta:
                    delta = minv[column]
                    column1 = column
            for column in range(n + 1):
                if used[column]:
                    u[p[column]] += delta
                    v[column] -= delta
                else:
                    minv[column] -= delta
            column0 = column1
            if p[column0] == 0:
                break
        while True:
            column1 = way[column0]
            p[column0] = p[column1]
            column0 = column1
            if column0 == 0:
                break
    answer = [-1] * n
    for column in range(1, n + 1):
        if p[column]:
            answer[p[column] - 1] = column - 1
    return answer


def support_indices(observed_logs, base_values, count):
    pairs = []
    for index, (observed, base) in enumerate(zip(observed_logs, base_values)):
        ratio = observed / base if base > 0.0 else math.inf
        pairs.append((index, ratio))
    pairs.sort(key=lambda item: (item[1], item[0]))
    best = None
    for start in range(0, len(pairs) - count + 1):
        window = pairs[start : start + count]
        logs = [math.log(item[1]) for item in window if item[1] > 0.0]
        if len(logs) != count:
            spread = math.inf
        else:
            spread = max(logs) - min(logs)
        candidate = (spread, tuple(item[0] for item in window))
        if best is None or candidate < best:
            best = candidate
    return best


class Reconciler:
    def __init__(self, payload):
        self.payload = payload
        self.histories = payload["reaction_histories"]
        self.bundles = payload["sample_bundles"]
        self.policy = payload["measurement_policy"]
        all_times = sorted(
            {
                float(observation["elapsed_s"])
                for bundle in self.bundles
                for observation in bundle["observations"]
            }
        )
        self.times = all_times
        self.time_index = {value: index for index, value in enumerate(all_times)}
        self.bundle_logs = [
            [
                -math.log1p(-float(observation["conversion"]))
                for observation in bundle["observations"]
            ]
            for bundle in self.bundles
        ]
        self.bundle_time_indices = [
            [
                self.time_index[float(observation["elapsed_s"])]
                for observation in bundle["observations"]
            ]
            for bundle in self.bundles
        ]

    def base_matrix(self, activation_energy):
        matrix = []
        for history in self.histories:
            scale = float(history["catalyst_scale"])
            matrix.append(
                [
                    scale * thermal_dose(history["steps"], elapsed, activation_energy)
                    for elapsed in self.times
                ]
            )
        return matrix

    def pair_evidence(self, activation_energy):
        bases = self.base_matrix(activation_energy)
        costs = []
        evidence = {}
        for history_index, history_base in enumerate(bases):
            row = []
            for bundle_index, logs in enumerate(self.bundle_logs):
                indices = self.bundle_time_indices[bundle_index]
                base_values = [history_base[index] for index in indices]
                count = len(logs) - 1
                usable = all(
                    value > float(self.policy["conversion_floor"])
                    for value in logs
                ) and all(value > 0.0 for value in base_values)
                if usable:
                    spread, support = support_indices(logs, base_values, count)
                else:
                    spread, support = 1e6, tuple(range(count))
                row.append(spread)
                evidence[(history_index, bundle_index)] = (spread, support)
            costs.append(row)
        return costs, evidence

    def assign(self, activation_energy):
        raw, evidence = self.pair_evidence(activation_energy)
        history_count = len(self.histories)
        bundle_count = len(self.bundles)
        size = max(history_count, bundle_count)
        dummy_cost = float(self.policy["conversion_tolerance"]) * 4.0
        costs = [[dummy_cost for _ in range(size)] for _ in range(size)]
        for history_index in range(history_count):
            for bundle_index in range(bundle_count):
                costs[history_index][bundle_index] = raw[history_index][bundle_index]
        for row in range(history_count, size):
            for column in range(bundle_count, size):
                costs[row][column] = 0.0
        assignment = hungarian(costs)
        assigned = {}
        for history_index in range(history_count):
            bundle_index = assignment[history_index]
            if bundle_index < bundle_count:
                assigned[history_index] = bundle_index
        objective = 0.0
        for history_index, bundle_index in assigned.items():
            objective += raw[history_index][bundle_index]
        objective += (history_count - len(assigned)) * dummy_cost
        return objective, assigned, evidence

    def estimate_activation_energy(self):
        low, high = (
            float(value)
            for value in self.policy["activation_energy_bounds_j_mol"]
        )
        sample_values = [
            low + (high - low) * index / 16.0
            for index in range(17)
        ]
        sampled = [
            (self.assign(value)[0], value)
            for value in sample_values
        ]
        best_index = min(range(len(sampled)), key=lambda index: (sampled[index][0], index))
        left = sample_values[max(0, best_index - 1)]
        right = sample_values[min(len(sample_values) - 1, best_index + 1)]
        ratio = (math.sqrt(5.0) - 1.0) / 2.0
        x1 = right - ratio * (right - left)
        x2 = left + ratio * (right - left)
        f1 = self.assign(x1)[0]
        f2 = self.assign(x2)[0]
        for _ in range(54):
            if f1 <= f2:
                right = x2
                x2, f2 = x1, f1
                x1 = right - ratio * (right - left)
                f1 = self.assign(x1)[0]
            else:
                left = x1
                x1, f1 = x2, f2
                x2 = left + ratio * (right - left)
                f2 = self.assign(x2)[0]
        candidates = [
            (sampled[best_index][0], sampled[best_index][1]),
            (f1, x1),
            (f2, x2),
        ]
        return min(candidates, key=lambda item: (item[0], item[1]))

    def solve(self):
        floor = float(self.policy["conversion_floor"])
        if any(
            max(logs) <= floor
            for logs in self.bundle_logs
        ):
            raise Refusal("no_consensus")
        _, activation_energy = self.estimate_activation_energy()
        objective, assigned, evidence = self.assign(activation_energy)
        if len(assigned) != len(self.histories):
            raise Refusal("no_consensus")

        inverse_values = []
        support_by_pair = {}
        for history_index, bundle_index in assigned.items():
            _, support = evidence[(history_index, bundle_index)]
            support_by_pair[(history_index, bundle_index)] = support
            base_indices = self.bundle_time_indices[bundle_index]
            history_base = self.base_matrix(activation_energy)[history_index]
            for point_index in support:
                base_value = history_base[base_indices[point_index]]
                inverse_values.append(
                    math.log(self.bundle_logs[bundle_index][point_index] / base_value)
                )
        inverse_values.sort()
        middle = len(inverse_values) // 2
        if len(inverse_values) % 2:
            ln_prefactor = inverse_values[middle]
        else:
            ln_prefactor = 0.5 * (
                inverse_values[middle - 1] + inverse_values[middle]
            )
        ln_low, ln_high = (
            float(value) for value in self.policy["ln_prefactor_bounds"]
        )
        if not ln_low <= ln_prefactor <= ln_high:
            raise Refusal("no_consensus")
        prefactor = math.exp(ln_prefactor)

        histories_by_index = {
            index: self.histories[index]["history_id"]
            for index in range(len(self.histories))
        }
        bundles_by_index = {
            index: self.bundles[index]["bundle_id"]
            for index in range(len(self.bundles))
        }
        assignments = []
        flagged = []
        residuals = []
        used_bundles = set()
        for history_index, bundle_index in sorted(
            assigned.items(), key=lambda item: histories_by_index[item[0]]
        ):
            used_bundles.add(bundle_index)
            support = list(support_by_pair[(history_index, bundle_index)])
            outliers = [
                index
                for index in range(len(self.bundles[bundle_index]["observations"]))
                if index not in support
            ]
            history_base = self.base_matrix(activation_energy)[history_index]
            for index in support:
                base_index = self.bundle_time_indices[bundle_index][index]
                predicted_log = prefactor * history_base[base_index]
                predicted = 1.0 - math.exp(-predicted_log)
                observed = self.bundles[bundle_index]["observations"][index]["conversion"]
                residual = predicted - float(observed)
                residuals.append(residual)
                if abs(residual) > float(self.policy["conversion_tolerance"]):
                    raise Refusal("no_consensus")
            for index in outliers:
                flagged.append(
                    {
                        "bundle_id": bundles_by_index[bundle_index],
                        "index": index,
                        "reason": "consensus_outlier",
                    }
                )
            assignments.append(
                {
                    "history_id": histories_by_index[history_index],
                    "bundle_id": bundles_by_index[bundle_index],
                    "supporting_indices": support,
                    "outlier_indices": outliers,
                }
            )
        for bundle_index in range(len(self.bundles)):
            if bundle_index in used_bundles:
                continue
            for index in range(len(self.bundles[bundle_index]["observations"])):
                flagged.append(
                    {
                        "bundle_id": bundles_by_index[bundle_index],
                        "index": index,
                        "reason": "no_consensus_history",
                    }
                )
        flagged.sort(key=lambda item: (item["bundle_id"], item["index"]))
        rmse = math.sqrt(
            math.fsum(residual * residual for residual in residuals) / len(residuals)
        )
        assignment_digest = hashlib.sha256(
            json.dumps(
                assignments,
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
            ).encode("utf-8")
        ).hexdigest()
        canonical_input = json.dumps(
            self.payload,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
        entry_id = hashlib.sha256(
            canonical_input
            + assignment_digest.encode("ascii")
            + f"{activation_energy:.6f}|{ln_prefactor:.6f}".encode("ascii")
        ).hexdigest()
        return {
            "schema_version": 1,
            "kinetics": {
                "activation_energy_j_mol": round(activation_energy, 6),
                "ln_prefactor": round(ln_prefactor, 6),
                "conversion_rmse": round(rmse, 8),
            },
            "assignments": assignments,
            "flagged_observations": flagged,
            "diagnostics": {
                "history_count": len(self.histories),
                "bundle_count": len(self.bundles),
                "assigned_count": len(assignments),
                "consensus_count": len(residuals),
                "unassigned_bundle_count": len(self.bundles) - len(assignments),
            },
            "ledger_entry_id": entry_id,
        }


def read_ledger(path):
    ledger_path = Path(path)
    if not ledger_path.exists():
        return []
    try:
        raw = ledger_path.read_bytes()
    except OSError as exc:
        raise Refusal("publish_failed") from exc
    if not raw:
        return []
    if not raw.endswith(b"\n"):
        raise Refusal("schema_invalid")
    entries = []
    for line in raw.splitlines():
        if not line or line.startswith(b" ") or line.endswith(b" "):
            raise Refusal("schema_invalid")
        try:
            entry = json.loads(line.decode("utf-8"), object_pairs_hook=_pairs)
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise Refusal("schema_invalid") from exc
        require_object(entry, LEDGER_KEYS)
        if require_int(entry["schema_version"]) != 1:
            raise Refusal("schema_invalid")
        for key in ("entry_id", "assignments_sha256"):
            value = require_string(entry[key])
            if len(value) != 64 or any(character not in "0123456789abcdef" for character in value):
                raise Refusal("schema_invalid")
        require_number(entry["activation_energy_j_mol"])
        require_number(entry["ln_prefactor"])
        entries.append(entry)
    return entries


def staged_write(path, payload):
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        dir=target.parent,
        prefix=f".{target.name}.",
    )
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        return Path(temporary)
    except Exception:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def publish(report, report_path, ledger_path):
    report_path = Path(report_path)
    ledger_path = Path(ledger_path)
    entries = read_ledger(ledger_path)
    entry = {
        "schema_version": 1,
        "entry_id": report["ledger_entry_id"],
        "activation_energy_j_mol": report["kinetics"]["activation_energy_j_mol"],
        "ln_prefactor": report["kinetics"]["ln_prefactor"],
        "assignments_sha256": hashlib.sha256(
            json.dumps(
                report["assignments"],
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
            ).encode("utf-8")
        ).hexdigest(),
    }
    if not entries or entries[-1]["entry_id"] != entry["entry_id"]:
        entries.append(entry)
    ledger_bytes = b"".join(
        json.dumps(item, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
        + b"\n"
        for item in entries
    )
    report_bytes = (
        json.dumps(report, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
        + b"\n"
    )

    try:
        old_report = report_path.read_bytes() if report_path.exists() else None
        old_ledger = ledger_path.read_bytes() if ledger_path.exists() else None
    except OSError as exc:
        raise Refusal("publish_failed") from exc
    report_temp = None
    ledger_temp = None
    report_replaced = False
    try:
        report_temp = staged_write(report_path, report_bytes)
        ledger_temp = staged_write(ledger_path, ledger_bytes)
        os.replace(report_temp, report_path)
        report_temp = None
        report_replaced = True
        os.replace(ledger_temp, ledger_path)
        ledger_temp = None
    except Exception as exc:
        if report_replaced:
            if old_report is None:
                try:
                    report_path.unlink()
                except OSError:
                    pass
            else:
                restore = staged_write(report_path, old_report)
                os.replace(restore, report_path)
        if old_ledger is not None and (
            not ledger_path.exists() or ledger_path.read_bytes() != old_ledger
        ):
            restore = staged_write(ledger_path, old_ledger)
            os.replace(restore, ledger_path)
        for temporary in (report_temp, ledger_temp):
            if temporary is not None:
                try:
                    temporary.unlink()
                except OSError:
                    pass
        raise Refusal("publish_failed") from exc


def run(input_path, output_path, ledger_path):
    payload = validate_dataset(read_json(input_path))
    report = Reconciler(payload).solve()
    publish(report, output_path, ledger_path)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--ledger", required=True)
    args = parser.parse_args(argv)
    try:
        run(args.input, args.output, args.ledger)
    except Refusal as exc:
        print(f"reconcile refused: {exc.code}", file=sys.stderr)
        return 3 if exc.code == "publish_failed" else 2
    except Exception:
        print("reconcile refused: internal_failure", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
