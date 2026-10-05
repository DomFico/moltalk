#!/usr/bin/env bash
# Write the tunnel-client profile: Secure MCP Tunnel -> moltalk over stdio (no network listener).
# Usage: scripts/configure-tunnel.sh tunnel_<32 hex chars>
set -euo pipefail
source "$(dirname "$0")/common.sh"
tunnel_id="${1:-}"
[[ "$tunnel_id" =~ ^tunnel_[0-9a-f]{32}$ ]] || { echo "Usage: $0 tunnel_<32 lowercase hex chars>" >&2; exit 2; }
[[ -x "$PROJECT_DIR/.venv/bin/moltalk" ]] || { echo "Missing $PROJECT_DIR/.venv; install the project first (see README)." >&2; exit 1; }
[[ -f "$KEY_FILE" ]] || { echo "Missing $KEY_FILE; run scripts/store-tunnel-key.sh first." >&2; exit 1; }
mkdir -p "$PROFILE_DIR"
"$TUNNEL_CLIENT_BIN" init --force --profile-dir "$PROFILE_DIR" --profile "$PROFILE" \
  --sample sample_mcp_stdio_local --tunnel-id "$tunnel_id" \
  --mcp-command "$PROJECT_DIR/.venv/bin/moltalk" \
  --control-plane-api-key-ref "file:$KEY_FILE" --health-listen-addr "$HEALTH_ADDR" >/dev/null
echo "Profile written to $PROFILE_DIR/$PROFILE.yaml"
"$TUNNEL_CLIENT_BIN" doctor --profile-dir "$PROFILE_DIR" --profile "$PROFILE"
