#!/usr/bin/env bash
# Remove the systemd user unit. Keeps the key and profile; delete ~/.config/moltalk to remove those.
set -euo pipefail
source "$(dirname "$0")/common.sh"
systemctl --user disable --now "$SERVICE" 2>/dev/null || true
rm -f "$HOME/.config/systemd/user/$SERVICE"
systemctl --user daemon-reload
echo "Removed $SERVICE."
