#!/usr/bin/env bash
# Show service state and tunnel-client readiness (local admin UI is loopback-only).
source "$(dirname "$0")/common.sh"
echo "service: $(systemctl --user is-active "$SERVICE" 2>/dev/null || echo inactive)"
printf 'healthz: '; curl -fsS "http://$HEALTH_ADDR/healthz" 2>/dev/null || echo "unreachable"; echo
printf 'readyz:  '; curl -fsS "http://$HEALTH_ADDR/readyz" 2>/dev/null || echo "not ready"; echo
echo "admin UI: http://$HEALTH_ADDR/ui    logs: journalctl --user -u $SERVICE -f"
