#!/usr/bin/env bash
# Stop the background tunnel. ChatGPT tool calls fail until it is started again.
set -euo pipefail
source "$(dirname "$0")/common.sh"
systemctl --user stop "$SERVICE" 2>/dev/null || true
echo "Stopped $SERVICE."
