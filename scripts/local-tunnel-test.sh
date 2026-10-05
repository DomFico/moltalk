#!/usr/bin/env bash
# Credential-free check of the exact tunnel-client -> stdio path, using tunnel-client's
# local in-memory control plane. Prints the local MCP URL for an MCP client to use.
set -euo pipefail
source "$(dirname "$0")/common.sh"
exec "$TUNNEL_CLIENT_BIN" dev proxy --mcp-command "$PROJECT_DIR/.venv/bin/moltalk" --print-json "$@"
