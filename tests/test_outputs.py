"""Behavioral acceptance for the polymerization reconciliation program."""

from __future__ import annotations

import copy
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pytest


APP = Path(os.environ.get("TASK_APP", "/app/polymer_reconcile.py"))
DATA = Path(os.environ.get("TASK_DATA_DIR", "/data"))
BASE_INPUT = DATA / "polymer_batches.json"
CALL_TIMEOUT_SEC = 45
PUBLIC_LIMIT_SEC = 40

EXPECTED_MAPPING = {
    "RX-11E-264B": "SPL-E09-2A79",
    "RX-14E-A610": "SPL-00E-C6A2",
    "RX-2D3-172D": "SPL-E3D-E560",
    "RX-3DA-7130": "SPL-8CC-A4D3",
    "RX-4B6-7C6B": "SPL-18B-0216",
    "RX-529-9A1A": "SPL-024-09A9",
    "RX-6E0-F554": "SPL-9D1-00CD",
    "RX-6F6-7B51": "SPL-20B-34D7",
    "RX-7B0-2578": "SPL-681-F63B",
    "RX-7E3-8B3F": "SPL-C12-A350",
    "RX-87B-7545": "SPL-01B-BA19",
    "RX-896-747A": "SPL-103-FD8C",
    "RX-997-5430": "SPL-2C6-3CEA",
    "RX-A17-514C": "SPL-B50-EBEF",
    "RX-A20-FB7F": "SPL-698-8754",
    "RX-B6F-DF3B": "SPL-094-0954",
    "RX-B72-79D9": "SPL-6CC-0D1E",
    "RX-C34-4046": "SPL-8FE-2484",
    "RX-C47-96B8": "SPL-692-A3AC",
    "RX-CE9-4CDD": "SPL-2F3-7801",
    "RX-D39-B38C": "SPL-B37-7BDD",
    "RX-E19-7357": "SPL-226-228D",
    "RX-E1C-EDAF": "SPL-2C7-F039",
    "RX-E8E-E27C": "SPL-95D-DD75",
    "RX-EBC-D174": "SPL-9E6-A5BD",
    "RX-F5B-C429": "SPL-E64-2CA0",
}
EXPECTED_OUTLIERS = {
    "RX-11E-264B": [0],
    "RX-14E-A610": [2],
    "RX-2D3-172D": [0],
    "RX-3DA-7130": [1],
    "RX-4B6-7C6B": [3],
    "RX-529-9A1A": [2],
    "RX-6E0-F554": [0],
    "RX-6F6-7B51": [1],
    "RX-7B0-2578": [3],
    "RX-7E3-8B3F": [0],
    "RX-87B-7545": [0],
    "RX-896-747A": [0],
    "RX-997-5430": [3],
    "RX-A17-514C": [0],
    "RX-A20-FB7F": [3],
    "RX-B6F-DF3B": [2],
    "RX-B72-79D9": [3],
    "RX-C34-4046": [3],
    "RX-C47-96B8": [2],
    "RX-CE9-4CDD": [1],
    "RX-D39-B38C": [1],
    "RX-E19-7357": [3],
    "RX-E1C-EDAF": [1],
    "RX-E8E-E27C": [2],
    "RX-EBC-D174": [1],
    "RX-F5B-C429": [2],
}
EXPECTED_UNASSIGNED = {"SPL-676-DB6F", "SPL-F21-6A18"}


def _invoke(
    input_path,
    output_path,
    ledger_path,
    timeout=CALL_TIMEOUT_SEC,
    app=APP,
):
    return subprocess.run(
        [
            sys.executable,
            str(app),
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--ledger",
            str(ledger_path),
        ],
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def _write_json(path, payload):
    path.write_text(
        json.dumps(payload, ensure_ascii=False, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )


def _base_payload():
    return json.loads(BASE_INPUT.read_text(encoding="utf-8"))


def _prepare_case(root, payload):
    input_path = root / "polymer_batches.json"
    output_path = root / "reconciliation.json"
    ledger_path = root / "experiment-ledger.jsonl"
    _write_json(input_path, payload)
    return input_path, output_path, ledger_path


def _subset_payload(count):
    payload = _base_payload()
    keep_histories = payload["reaction_histories"][:count]
    keep_ids = {item["history_id"] for item in keep_histories}
    keep_bundles = [
        item
        for item in payload["sample_bundles"]
        if item["bundle_id"] in {EXPECTED_MAPPING[history_id] for history_id in keep_ids}
    ]
    payload["reaction_histories"] = keep_histories
    payload["sample_bundles"] = keep_bundles
    payload["measurement_policy"]["max_assignments"] = count
    return payload


@pytest.fixture(scope="session")
def lab_case(tmp_path_factory):
    root = tmp_path_factory.mktemp("laboratory")
    input_path = root / "polymer_batches.json"
    output_path = root / "reconciliation.json"
    ledger_path = root / "experiment-ledger.jsonl"
    shutil.copyfile(BASE_INPUT, input_path)
    process = _invoke(input_path, output_path, ledger_path)
    return {
        "root": root,
        "input": input_path,
        "output": output_path,
        "ledger": ledger_path,
        "run": process,
    }


def test_report_file_is_published(lab_case):
    """题面：「报告文件必须是紧凑 JSON」。"""
    assert lab_case["output"].is_file()


def test_success_is_silent(lab_case):
    """题面：「成功时标准输出和标准错误都为空。」"""
    run = lab_case["run"]
    assert run.stdout == ""
    assert run.stderr == ""


def test_natural_gate_reconstructs_reaction_assignments(lab_case):
    """题面：「每个反应釜和托盘都必须在各自方向上一一对应」。"""
    report = json.loads(lab_case["output"].read_text(encoding="utf-8"))
    actual = {
        item["history_id"]: item["bundle_id"] for item in report["assignments"]
    }
    assert actual == EXPECTED_MAPPING


def test_report_has_exact_top_level_keys(lab_case):
    """题面：「顺序也按这里列出的顺序写出。」"""
    report = json.loads(lab_case["output"].read_text(encoding="utf-8"))
    assert set(report) == {
        "schema_version",
        "kinetics",
        "assignments",
        "flagged_observations",
        "diagnostics",
        "ledger_entry_id",
    }


def test_report_preserves_declared_key_order(lab_case):
    """题面：「顺序也按这里列出的顺序写出。」"""
    report = json.loads(lab_case["output"].read_text(encoding="utf-8"))
    assert list(report) == [
        "schema_version",
        "kinetics",
        "assignments",
        "flagged_observations",
        "diagnostics",
        "ledger_entry_id",
    ]


def test_kinetics_has_exact_keys(lab_case):
    """[REPORT-02] [KIN-KEYS-02] [KIN-KEYS-03] 题面：「前两项分别采用焦耳每摩尔和对数形式，均为有限数值。」"""
    report = json.loads(lab_case["output"].read_text(encoding="utf-8"))
    assert set(report["kinetics"]) == {
        "activation_energy_j_mol",
        "ln_prefactor",
        "conversion_rmse",
    }


def test_assignment_items_have_exact_keys(lab_case):
    """题面：「`assignments` 是对象数组，按 `history_id` 的 UTF-8 字节序升序排列；每个对象恰好有 `history_id`、`bundle_id`、`supporting_indices`、`outlier_indices` 四个键。」"""
    report = json.loads(lab_case["output"].read_text(encoding="utf-8"))
    assert all(
        set(item)
        == {
            "history_id",
            "bundle_id",
            "supporting_indices",
            "outlier_indices",
        }
        for item in report["assignments"]
    )


def test_flagged_items_have_exact_keys(lab_case):
    """题面：「每项恰好有 `bundle_id`、`index`、`reason` 三个键。」"""
    report = json.loads(lab_case["output"].read_text(encoding="utf-8"))
    assert all(
        set(item) == {"bundle_id", "index", "reason"}
        for item in report["flagged_observations"]
    )


def test_diagnostics_has_exact_keys(lab_case):
    """题面：「`diagnostics` 恰好有 `history_count`、`bundle_count`、`assigned_count`、`consensus_count`、`unassigned_bundle_count` 五个键」。"""
    report = json.loads(lab_case["output"].read_text(encoding="utf-8"))
    assert set(report["diagnostics"]) == {
        "history_count",
        "bundle_count",
        "assigned_count",
        "consensus_count",
        "unassigned_bundle_count",
    }


def test_flagged_observations_have_declared_order_and_reasons(lab_case):
    """题面：「`flagged_observations` 按 `bundle_id` 的 UTF-8 字节序升序，同一托盘再按 `index` 升序。」"""
    report = json.loads(lab_case["output"].read_text(encoding="utf-8"))
    actual = [
        (item["bundle_id"], item["index"])
        for item in report["flagged_observations"]
    ]
    assert actual == sorted(actual)
    reasons = {
        item["reason"] for item in report["flagged_observations"]
    }
    assert reasons == {"consensus_outlier", "no_consensus_history"}


def test_ledger_has_exact_line_keys(lab_case):
    """[LEDGER-01] [LEDGER-06] 题面：「每一行恰好有 `schema_version`、`entry_id`、`activation_energy_j_mol`、`ln_prefactor`、`assignments_sha256` 五个键」。"""
    lines = lab_case["ledger"].read_text(encoding="utf-8").splitlines()
    entry = json.loads(lines[0])
    assert list(entry) == [
        "schema_version",
        "entry_id",
        "activation_energy_j_mol",
        "ln_prefactor",
        "assignments_sha256",
    ]


def test_report_is_compact_json_with_one_lf(lab_case):
    """题面：「报告文件必须是紧凑 JSON：不使用缩进或分隔空格」。"""
    raw = lab_case["output"].read_bytes()
    assert raw.endswith(b"\n")
    assert not raw.endswith(b"\n\n")
    assert b": " not in raw
    assert b", " not in raw


def test_ledger_is_compact_jsonl_with_one_lf_per_entry(lab_case):
    """题面：「行内用紧凑 JSON，每行一个 LF，键顺序与上面相同。」"""
    raw = lab_case["ledger"].read_bytes()
    assert raw.endswith(b"\n")
    assert len(raw.splitlines()) == 1
    assert b": " not in raw
    assert b", " not in raw


def test_report_uses_six_and_eight_decimal_kinetics(lab_case):
    """[REPORT-07] [REPORT-TIE-01] 题面：「报告数值按十进制四舍五入保留六位小数」。"""
    raw = lab_case["output"].read_text(encoding="utf-8")
    assert re.search(r'"activation_energy_j_mol":\d+\.\d{6}', raw)
    assert re.search(r'"ln_prefactor":\d+\.\d{6}', raw)
    assert re.search(r'"conversion_rmse":\d+\.\d{8}', raw)


def test_assignments_are_sorted_by_history_id(lab_case):
    """题面：「`assignments` 是对象数组，按 `history_id` 的 UTF-8 字节序升序排列；」。"""
    report = json.loads(lab_case["output"].read_text(encoding="utf-8"))
    actual = [item["history_id"] for item in report["assignments"]]
    assert actual == sorted(actual, key=lambda value: value.encode("utf-8"))


def test_every_assignment_covers_all_observation_indices(lab_case):
    """[REPORT-09] [ASSIGN-KEYS-02] 题面：「两个索引数组都是非负整数数组、严格升序，且合起来恰好覆盖该托盘的全部读数下标」。"""
    report = json.loads(lab_case["output"].read_text(encoding="utf-8"))
    assert all(
        sorted(item["supporting_indices"] + item["outlier_indices"]) == [0, 1, 2, 3]
        for item in report["assignments"]
    )


def test_assignments_match_laboratory_mapping(lab_case):
    """题面：「每个反应釜和托盘都必须在各自方向上一一对应」。"""
    report = json.loads(lab_case["output"].read_text(encoding="utf-8"))
    actual = {
        item["history_id"]: item["bundle_id"] for item in report["assignments"]
    }
    assert actual == EXPECTED_MAPPING


def test_outlier_consensus_matches_laboratory_notes(lab_case):
    """题面：「每个对应至少保留三个支持点，最多有一个异常点。」"""
    report = json.loads(lab_case["output"].read_text(encoding="utf-8"))
    actual = {
        item["history_id"]: item["outlier_indices"] for item in report["assignments"]
    }
    assert actual == EXPECTED_OUTLIERS


def test_unassigned_bundle_readings_are_flagged(lab_case):
    """题面：「无对应反应釜的托盘所有读数写 `no_consensus_history`。」"""
    report = json.loads(lab_case["output"].read_text(encoding="utf-8"))
    actual = {
        item["bundle_id"]
        for item in report["flagged_observations"]
        if item["reason"] == "no_consensus_history"
    }
    assert actual == EXPECTED_UNASSIGNED


def test_kinetics_match_laboratory_reference(lab_case):
    """题面：「复核结果要求活化能误差不超过 600 J/mol，对数指前因子误差不超过 0.020」。"""
    report = json.loads(lab_case["output"].read_text(encoding="utf-8"))
    assert abs(report["kinetics"]["activation_energy_j_mol"] - 92500.0) <= 600.0
    assert abs(report["kinetics"]["ln_prefactor"] - 20.888673) <= 0.020
    assert report["kinetics"]["conversion_rmse"] <= 0.0015


def test_diagnostics_count_the_laboratory_fixture(lab_case):
    """[REPORT-10] [DIAG-COUNT-01] [DIAG-02] 题面：「`assigned_count` 必须等于反应釜数量，`unassigned_bundle_count` 必须等于托盘总数减已对应托盘数」。"""
    report = json.loads(lab_case["output"].read_text(encoding="utf-8"))
    assert report["diagnostics"] == {
        "history_count": 26,
        "bundle_count": 28,
        "assigned_count": 26,
        "consensus_count": 78,
        "unassigned_bundle_count": 2,
    }


def test_ledger_entry_id_is_lowercase_sha256(lab_case):
    """[REPORT-11] [LEDGER-07] 题面：「`ledger_entry_id` 是 64 位小写十六进制字符串。」"""
    report = json.loads(lab_case["output"].read_text(encoding="utf-8"))
    assert re.fullmatch(r"[0-9a-f]{64}", report["ledger_entry_id"])


def test_repeated_input_keeps_report_bytes(lab_case):
    """题面：「连续处理同一份输入两次，报告和交接记录都必须逐字节不变。」"""
    before = lab_case["output"].read_bytes()
    process = _invoke(lab_case["input"], lab_case["output"], lab_case["ledger"])
    assert process.returncode == 0
    assert lab_case["output"].read_bytes() == before


def test_repeated_input_keeps_ledger_bytes(lab_case):
    """题面：「重复同一输入时不追加。」"""
    before = lab_case["ledger"].read_bytes()
    process = _invoke(lab_case["input"], lab_case["output"], lab_case["ledger"])
    assert process.returncode == 0
    assert lab_case["ledger"].read_bytes() == before


def test_missing_input_refuses_without_output(tmp_path):
    """题面：「输入不存在、无法读取或 UTF-8 解码失败时，退出码为 2」。"""
    output = tmp_path / "reconciliation.json"
    ledger = tmp_path / "experiment-ledger.jsonl"
    process = _invoke(tmp_path / "not-present.json", output, ledger)
    assert process.returncode == 2
    assert process.stdout == ""
    assert process.stderr == "reconcile refused: input_unavailable\n"
    assert not output.exists()
    assert not ledger.exists()


def test_json_parse_failure_leaves_targets_absent(tmp_path):
    """题面：「输入不存在、无法读取或 UTF-8 解码失败时，退出码为 2，标准输出为空，标准错误只写一行 `reconcile refused: input_unavailable` 或 `reconcile refused: json_invalid`。」"""
    source = tmp_path / "payload_case.json"
    source.write_text("{broken", encoding="utf-8")
    output = tmp_path / "reconciliation.json"
    ledger = tmp_path / "experiment-ledger.jsonl"
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: json_invalid\n"
    assert not output.exists()
    assert not ledger.exists()


def test_duplicate_top_level_key_is_refused(tmp_path):
    """[JSON-05] [TOP-SET-02] 题面：「根对象除 schema_version、reaction_histories、sample_bundles、measurement_policy 外不得再放成员。」；「同一对象内重复键都必须拒绝。」"""
    raw = BASE_INPUT.read_text(encoding="utf-8")
    duplicate = raw.replace(
        '{\n  "schema_version": 1,',
        '{\n  "schema_version": 1,\n  "schema_version": 1,',
        1,
    )
    source = tmp_path / "duplicate.json"
    source.write_text(duplicate, encoding="utf-8")
    output = tmp_path / "reconciliation.json"
    ledger = tmp_path / "experiment-ledger.jsonl"
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: duplicate_key\n"
    assert not output.exists()


def test_duplicate_nested_key_is_refused(tmp_path):
    """题面：「同一对象内重复键都必须拒绝。」"""
    raw = BASE_INPUT.read_text(encoding="utf-8")
    duplicate = raw.replace(
        '"history_id": "RX-11E-264B"',
        '"history_id": "RX-11E-264B", "history_id": "RX-11E-264B"',
        1,
    )
    source = tmp_path / "nested.json"
    source.write_text(duplicate, encoding="utf-8")
    output = tmp_path / "reconciliation.json"
    ledger = tmp_path / "experiment-ledger.jsonl"
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: duplicate_key\n"
    assert not output.exists()


@pytest.mark.parametrize("mutation", ["missing", "extra"])
def test_top_level_key_set_is_exact(tmp_path, mutation):
    """题面：「根对象除 schema_version、reaction_histories、sample_bundles、measurement_policy 外不得再放成员。」"""
    payload = _base_payload()
    mutation_action = {
        "missing": lambda value: value.pop("measurement_policy"),
        "extra": lambda value: value.__setitem__("unexpected", 1),
    }
    mutation_action[mutation](payload)
    source, output, ledger = _prepare_case(tmp_path, payload)
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: schema_invalid\n"
    assert not output.exists()


@pytest.mark.parametrize("field", ["history_id", "steps", "catalyst_scale"])
def test_history_key_set_is_exact(tmp_path, field):
    """题面：「每个对象恰好有 `history_id`、`steps`、`catalyst_scale`。」"""
    payload = _base_payload()
    payload["reaction_histories"][0].pop(field)
    source, output, ledger = _prepare_case(tmp_path, payload)
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: schema_invalid\n"
    assert not output.exists()


@pytest.mark.parametrize("field", ["bundle_id", "observations"])
def test_bundle_key_set_is_exact(tmp_path, field):
    """题面：「每个对象的键恰好是 `bundle_id`、`observations`。」"""
    payload = _base_payload()
    payload["sample_bundles"][0].pop(field)
    source, output, ledger = _prepare_case(tmp_path, payload)
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: schema_invalid\n"
    assert not output.exists()


@pytest.mark.parametrize("field", ["elapsed_s", "conversion"])
def test_observation_key_set_is_exact(tmp_path, field):
    """题面：「每个对象恰好有 `elapsed_s`、`conversion`；」。"""
    payload = _base_payload()
    payload["sample_bundles"][0]["observations"][0].pop(field)
    source, output, ledger = _prepare_case(tmp_path, payload)
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: schema_invalid\n"
    assert not output.exists()


def test_boolean_is_not_a_numeric_measurement(tmp_path):
    """[TYPE-01] [ERROR-01] 题面：「布尔值和数字字符串都不算数值。」"""
    payload = _base_payload()
    payload["reaction_histories"][0]["catalyst_scale"] = True
    source, output, ledger = _prepare_case(tmp_path, payload)
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: schema_invalid\n"
    assert not output.exists()


def test_schema_version_boolean_is_refused(tmp_path):
    """题面：「`schema_version` 必须是非布尔整数 `1`。」"""
    payload = _base_payload()
    payload["schema_version"] = True
    source, output, ledger = _prepare_case(tmp_path, payload)
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: schema_invalid\n"
    assert not output.exists()


def test_duplicate_history_identifier_is_refused(tmp_path):
    """[HISTORY-01] [HIST-OBJ-02] 题面：「`history_id` 是非空字符串，去首尾空白后长度 6 到 24，批次内不能重复。」"""
    payload = _base_payload()
    payload["reaction_histories"][1]["history_id"] = payload["reaction_histories"][0]["history_id"]
    source, output, ledger = _prepare_case(tmp_path, payload)
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: constraint_invalid\n"
    assert not output.exists()


def test_catalyst_scale_above_public_bound_is_refused(tmp_path):
    """题面：「`catalyst_scale` 是大于等于 0.2 且小于等于 3.0 的有限数值。」"""
    payload = _base_payload()
    payload["reaction_histories"][0]["catalyst_scale"] = 3.01
    source, output, ledger = _prepare_case(tmp_path, payload)
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: constraint_invalid\n"
    assert not output.exists()


def test_step_object_has_exact_keys(tmp_path):
    """题面：「`steps` 是 2 到 8 个对象的数组，每个对象的键恰好是 `duration_s`、`temp_start_k`、`temp_end_k`。」"""
    payload = _base_payload()
    payload["reaction_histories"][0]["steps"][0]["ramp_rate"] = 1.0
    source, output, ledger = _prepare_case(tmp_path, payload)
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: schema_invalid\n"
    assert not output.exists()


@pytest.mark.parametrize("duration", [0, 7200.01])
def test_step_duration_range_is_enforced(tmp_path, duration):
    """题面：「`duration_s` 必须大于 0 且不超过 7200。」"""
    payload = _base_payload()
    payload["reaction_histories"][0]["steps"][0]["duration_s"] = duration
    source, output, ledger = _prepare_case(tmp_path, payload)
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: constraint_invalid\n"
    assert not output.exists()


def test_step_temperature_above_public_bound_is_refused(tmp_path):
    """题面：「两个温度是 280 到 500 K 的有限数值。」"""
    payload = _base_payload()
    payload["reaction_histories"][0]["steps"][0]["temp_end_k"] = 500.1
    source, output, ledger = _prepare_case(tmp_path, payload)
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: constraint_invalid\n"
    assert not output.exists()


def test_observation_count_below_public_minimum_is_refused(tmp_path):
    """[OBS-01] [STEP-ARR-02] 题面：「`observations` 是 4 到 6 个对象的数组。」"""
    payload = _base_payload()
    payload["sample_bundles"][0]["observations"].pop()
    source, output, ledger = _prepare_case(tmp_path, payload)
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: constraint_invalid\n"
    assert not output.exists()


def test_measurement_policy_has_exact_keys(tmp_path):
    """[POLICY-01] [POLICY-03] 题面：「`measurement_policy` 恰好有这六个键：」"""
    payload = _base_payload()
    payload["measurement_policy"]["unknown_rule"] = 1
    source, output, ledger = _prepare_case(tmp_path, payload)
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: schema_invalid\n"
    assert not output.exists()


def test_max_assignments_must_equal_history_count(tmp_path):
    """题面：「`max_assignments` | 非布尔整数，必须等于反应釜数量，且不超过 30」"""
    payload = _base_payload()
    payload["measurement_policy"]["max_assignments"] = 25
    source, output, ledger = _prepare_case(tmp_path, payload)
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: constraint_invalid\n"
    assert not output.exists()


def test_conversion_one_is_outside_the_public_range(tmp_path):
    """[RANGE-01] [ERROR-01] 题面：「`conversion` 是大于等于 0 且严格小于 1 的有限数。」"""
    payload = _base_payload()
    payload["sample_bundles"][0]["observations"][0]["conversion"] = 1.0
    source, output, ledger = _prepare_case(tmp_path, payload)
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: constraint_invalid\n"
    assert not output.exists()


def test_elapsed_zero_is_outside_the_public_range(tmp_path):
    """[RANGE-02] [ERROR-03] 题面：「`elapsed_s` 是大于 0 且不超过 10000 的有限秒数。」"""
    payload = _base_payload()
    payload["sample_bundles"][0]["observations"][0]["elapsed_s"] = 0
    source, output, ledger = _prepare_case(tmp_path, payload)
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: constraint_invalid\n"
    assert not output.exists()


def test_public_boundary_accepts_eight_histories(tmp_path):
    """题面：「`reaction_histories` 是 8 到 30 个对象的数组。」"""
    payload = _subset_payload(8)
    source, output, ledger = _prepare_case(tmp_path, payload)
    process = _invoke(source, output, ledger)
    assert process.returncode == 0
    assert output.exists()
    assert ledger.exists()


def test_public_boundary_rejects_seven_histories(tmp_path):
    """题面：「`reaction_histories` 是 8 到 30 个对象的数组。」"""
    payload = _subset_payload(8)
    payload["reaction_histories"].pop()
    payload["measurement_policy"]["max_assignments"] = 7
    source, output, ledger = _prepare_case(tmp_path, payload)
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: constraint_invalid\n"
    assert not output.exists()


def test_public_boundary_accepts_maximum_batch_within_limit(tmp_path):
    """题面：「整批最多 30 个反应釜、32 个样品托盘、每个托盘 6 个读数时，程序应在 40 秒内完成。」"""
    source = tmp_path / "polymer_batches.json"
    output = tmp_path / "reconciliation.json"
    ledger = tmp_path / "experiment-ledger.jsonl"
    shutil.copyfile(DATA / "polymer_batches_max.json", source)
    started = time.monotonic()
    process = _invoke(source, output, ledger, timeout=45)
    elapsed = time.monotonic() - started
    assert process.returncode == 0
    assert output.exists()
    assert elapsed < 40


def test_public_boundary_rejects_31_histories(tmp_path):
    """题面：「`reaction_histories` 是 8 到 30 个对象的数组。」"""
    payload = _base_payload()
    for offset in range(5):
        history = copy.deepcopy(payload["reaction_histories"][offset])
        history["history_id"] = f"RX-OVER-{offset:02d}"
        payload["reaction_histories"].append(history)
    payload["measurement_policy"]["max_assignments"] = 31
    source, output, ledger = _prepare_case(tmp_path, payload)
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: constraint_invalid\n"
    assert not output.exists()


def test_public_boundary_rejects_33_bundles(tmp_path):
    """题面：「`sample_bundles` 是 8 到 32 个对象的数组。」"""
    payload = _base_payload()
    for offset in range(5):
        bundle = copy.deepcopy(payload["sample_bundles"][offset])
        bundle["bundle_id"] = f"SPL-OVER-{offset:02d}"
        payload["sample_bundles"].append(bundle)
    source, output, ledger = _prepare_case(tmp_path, payload)
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: constraint_invalid\n"
    assert not output.exists()


def test_changed_input_appends_one_ledger_entry(tmp_path):
    """题面：「输入内容改变并产生不同复核结果时，追加一行」。"""
    payload = _subset_payload(8)
    source, output, ledger = _prepare_case(tmp_path, payload)
    first = _invoke(source, output, ledger)
    payload["reaction_histories"][0]["steps"][0]["duration_s"] += 25.0
    _write_json(source, payload)
    second = _invoke(source, output, ledger)
    assert first.returncode == 0
    assert second.returncode == 0
    assert len(ledger.read_text(encoding="utf-8").splitlines()) == 2


def test_invalid_existing_ledger_preserves_both_files(tmp_path):
    """题面：「已有报告和交接记录的字节必须保持不变」。"""
    payload = _subset_payload(8)
    source, output, ledger = _prepare_case(tmp_path, payload)
    output.write_bytes(b"old-report\n")
    ledger.write_bytes(b"{}\n")
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: schema_invalid\n"
    assert output.read_bytes() == b"old-report\n"
    assert ledger.read_bytes() == b"{}\n"


def test_failure_does_not_rewrite_the_input(tmp_path):
    """题面：「已有报告和交接记录的字节必须保持不变，输入文件也不得被改写。」"""
    payload = _subset_payload(8)
    payload["sample_bundles"][0]["observations"][0]["conversion"] = 1.0
    source, output, ledger = _prepare_case(tmp_path, payload)
    before = source.read_bytes()
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert source.read_bytes() == before
    assert not output.exists()


def test_failed_publish_restores_existing_ledger_bytes(tmp_path):
    """[FAIL-02] [ERROR-02] 题面：「发布失败可能发生在报告已经换过一次之后，因此恢复时要把两个旧文件都放回原字节。」"""
    payload = _subset_payload(8)
    source, output, ledger = _prepare_case(tmp_path, payload)
    first = _invoke(source, output, ledger)
    old_report = output.read_bytes()
    old_ledger = ledger.read_bytes()
    payload["reaction_histories"][0]["steps"][0]["duration_s"] += 20.0
    _write_json(source, payload)
    wrapper = tmp_path / "second-replace-failure.py"
    wrapper.write_text(
        "\n".join(
            [
                "import importlib.util",
                "import os",
                "import sys",
                f"app = {str(APP)!r}",
                "spec = importlib.util.spec_from_file_location('candidate', app)",
                "candidate = importlib.util.module_from_spec(spec)",
                "spec.loader.exec_module(candidate)",
                "real_replace = os.replace",
                "replace_count = 0",
                "def fail_second_replace(source, target):",
                "    global replace_count",
                "    replace_count += 1",
                "    if replace_count == 2:",
                "        raise OSError('injected second replace failure')",
                "    return real_replace(source, target)",
                "candidate.os.replace = fail_second_replace",
                "raise SystemExit(candidate.main(sys.argv[1:]))",
                "",
            ]
        ),
        encoding="utf-8",
    )
    process = _invoke(source, output, ledger, app=wrapper)
    assert first.returncode == 0
    assert process.returncode == 3
    assert process.stderr == "reconcile refused: publish_failed\n"
    assert output.read_bytes() == old_report
    assert ledger.read_bytes() == old_ledger
    assert list(tmp_path.glob(".reconciliation.json.*")) == []
    assert list(tmp_path.glob(".experiment-ledger.jsonl.*")) == []


def test_inconsistent_bundle_refuses_with_no_consensus(tmp_path):
    """题面：「没有唯一可接受的一致对应时，退出码为 2」。"""
    payload = _subset_payload(8)
    for bundle in payload["sample_bundles"]:
        for observation in bundle["observations"]:
            observation["conversion"] = 0.0
    source, output, ledger = _prepare_case(tmp_path, payload)
    process = _invoke(source, output, ledger)
    assert process.returncode == 2
    assert process.stderr == "reconcile refused: no_consensus\n"
    assert not output.exists()
    assert not ledger.exists()
