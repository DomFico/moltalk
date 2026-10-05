# Shared paths for the tunnel scripts. Override any of these in the environment.
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# scripts/fetch-tools.sh puts tunnel-client in tools/; an older layout kept it next to the project.
if [[ -z "${TUNNEL_CLIENT_BIN:-}" ]]; then
  TUNNEL_CLIENT_BIN="$PROJECT_DIR/tools/tunnel-client/tunnel-client"
  [[ -x "$TUNNEL_CLIENT_BIN" ]] || TUNNEL_CLIENT_BIN="$PROJECT_DIR/../tools/tunnel-client/tunnel-client"
fi
CONFIG_DIR="${MOLTALK_CONFIG_DIR:-$HOME/.config/moltalk}"
KEY_FILE="$CONFIG_DIR/tunnel-runtime-key"
PROFILE_DIR="$CONFIG_DIR/tunnel-profiles"
PROFILE=moltalk
HEALTH_ADDR="${MOLTALK_HEALTH_ADDR:-127.0.0.1:8791}"
SERVICE=moltalk-tunnel.service
