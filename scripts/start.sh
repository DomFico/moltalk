#!/usr/bin/env bash
# Start the tunnel. Default: background systemd user service. Use --foreground to run in this terminal.
set -euo pipefail
source "$(dirname "$0")/common.sh"
[[ -f "$PROFILE_DIR/$PROFILE.yaml" ]] || { echo "No profile; run scripts/configure-tunnel.sh first." >&2; exit 1; }
if [[ "${1:-}" == "--foreground" ]]; then
  exec "$TUNNEL_CLIENT_BIN" run --profile-dir "$PROFILE_DIR" --profile "$PROFILE"
fi
unit_dir="$HOME/.config/systemd/user"
mkdir -p "$unit_dir"
cat > "$unit_dir/$SERVICE" <<UNIT
[Unit]
Description=MolTalk MCP server via OpenAI Secure MCP Tunnel
After=network-online.target

[Service]
ExecStart=$TUNNEL_CLIENT_BIN run --profile-dir $PROFILE_DIR --profile $PROFILE
Restart=on-failure
RestartSec=5
MemoryMax=3G
TasksMax=256
NoNewPrivileges=yes

[Install]
WantedBy=default.target
UNIT
systemctl --user daemon-reload
systemctl --user restart "$SERVICE"
for _ in $(seq 1 30); do
  curl -fsS "http://$HEALTH_ADDR/readyz" >/dev/null 2>&1 && break
  sleep 1
done
"$(dirname "$0")/status.sh"
