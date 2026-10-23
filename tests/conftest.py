"""Expose stable behavior case identifiers in the verifier transcript."""

from __future__ import annotations

import inspect

import pytest


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    outcome = yield
    report = outcome.get_result()
    if report.when == "call":
        item._behavior_case_outcome = report.outcome


def pytest_sessionfinish(session, exitstatus):
    for item in session.items:
        doc = inspect.getdoc(item.obj) or item.name
        outcome = getattr(item, "_behavior_case_outcome", "failed")
        print(f"OBSERVED_CASE {doc} {outcome.upper()}")
