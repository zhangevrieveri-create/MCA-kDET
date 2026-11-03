#!/usr/bin/env bash
# Runs an Rscript command with a watchdog against the intermittent R hangs
# (CLAUDE.md, "已知问题"): if the console log has not grown for STALL seconds
# (startup hang: STARTUP seconds), the run is killed and restarted; run.R skips
# the chunks already on disk. Gives up after TRIES attempts.
#
#   bash tests/run_with_watchdog.sh <console_log> <done_regex> -- <command...>
# Environment: STARTUP (default 120), STALL (default 900), TRIES (default 6).
log=$1; done_re=$2; shift 3
STARTUP=${STARTUP:-120}; STALL=${STALL:-900}; TRIES=${TRIES:-6}
for a in $(seq 1 "$TRIES"); do
  echo "=== attempt $a: $(date '+%F %T')" >> "$log.attempts"
  : > "$log"
  "$@" > "$log" 2>&1 &
  pid=$!; last=0; still=0; limit=$STARTUP
  while kill -0 "$pid" 2>/dev/null; do
    sleep 10
    sz=$(stat -c %s "$log" 2>/dev/null || echo 0)
    if [ "$sz" -gt "$last" ]; then last=$sz; still=0; limit=$STALL; else still=$((still + 10)); fi
    if [ "$still" -ge "$limit" ]; then
      echo "    no output for ${still}s, killing pid $pid" >> "$log.attempts"
      kill "$pid" 2>/dev/null; sleep 5; kill -9 "$pid" 2>/dev/null
      break
    fi
  done
  wait "$pid" 2>/dev/null; rc=$?
  echo "    exit $rc" >> "$log.attempts"
  if grep -qE "$done_re" "$log"; then echo "DONE after $a attempt(s)"; exit 0; fi
  sleep 5
done
echo "GAVE UP after $TRIES attempts"; exit 1
