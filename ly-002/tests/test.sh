#!/bin/bash
set -euo pipefail

LOG_DIR=/logs/verifier
# Fail-closed reward path: /logs/verifier/reward.txt.
REWARD_FILE="$LOG_DIR/reward.txt"
STATUS_FILE="$LOG_DIR/verification-status.txt"
TEST_LOG="$LOG_DIR/test-stdout.txt"
SCORING_LOG="$LOG_DIR/scoring-report.json"
REPORT="$LOG_DIR/pytest-results.xml"
completed=0

mkdir -p "$LOG_DIR"

write_reward() {
  chmod 0644 "$REWARD_FILE" 2>/dev/null || true
  printf '%s\n' "$1" > "$LOG_DIR/reward.txt"
  chmod 0444 "$REWARD_FILE"
}

write_status() {
  chmod 0644 "$STATUS_FILE" 2>/dev/null || true
  printf '%s\n' "$1" > "$STATUS_FILE"
  chmod 0444 "$STATUS_FILE"
}

on_exit() {
  if [ "$completed" -ne 1 ]; then
    write_status "INCOMPLETE: verifier exited before pytest completion"
    write_reward 0.000
  fi
}

on_signal() {
  write_status "INCOMPLETE: verifier interrupted before pytest completion"
  write_reward 0.000
  completed=1
  exit 0
}

write_reward 0.000
write_status "INCOMPLETE: pytest has not completed"
trap on_exit EXIT
trap on_signal INT TERM HUP

set +e
timeout --signal=TERM --kill-after=10s 60s \
  python3 -m pip install --no-index --no-cache-dir \
  --find-links=/tests/wheels pytest==8.4.1 > "$TEST_LOG" 2>&1
INSTALL_RC=$?
set -e
if [ "$INSTALL_RC" -ne 0 ]; then
  cat "$TEST_LOG" 2>/dev/null || true
  write_status "INCOMPLETE: local test dependency installation failed"
  write_reward 0.000
  completed=1
  exit 0
fi

rm -f "$REPORT"
set +e
env -u PYTEST_ADDOPTS -u PYTEST_PLUGINS \
  timeout --signal=TERM --kill-after=10s 300s \
  python3 -u -m pytest -v -s -p no:cacheprovider \
  --junitxml="$REPORT" /tests/test_outputs.py >> "$TEST_LOG" 2>&1
PYTEST_RC=$?
set -e
cat "$TEST_LOG" 2>/dev/null || true

if [ "$PYTEST_RC" -ne 0 ] && [ "$PYTEST_RC" -ne 1 ]; then
  case "$PYTEST_RC" in
    124|125|126|127|137|143)
      write_status "INCOMPLETE: pytest did not complete; runner exited with code $PYTEST_RC"
      ;;
    *)
      write_status "INCOMPLETE: pytest exited before a scoreable report with code $PYTEST_RC"
      ;;
  esac
  write_reward 0.000
elif [ ! -s "$REPORT" ]; then
  write_status "INCOMPLETE: pytest completion report missing"
  write_reward 0.000
else
  if python3 /tests/score.py \
    --config /tests/scoring.json \
    --junit "$REPORT" \
    --report-out "$SCORING_LOG" \
    --reward-out "$REWARD_FILE"; then
    REWARD=$(cat "$REWARD_FILE")
    case "$REWARD" in
      0.[0-9][0-9][0-9]|1.000) ;;
      *)
        write_status "INCOMPLETE: weighted scorer produced an invalid reward"
        write_reward 0.000
        completed=1
        exit 0
        ;;
    esac
    chmod 0444 "$REWARD_FILE"
    write_status "COMPLETED: weighted scoring finished with reward=$REWARD"
  else
    write_status "INCOMPLETE: weighted scoring rejected the report"
    write_reward 0.000
  fi
fi

completed=1
exit 0
