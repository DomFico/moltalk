#!/usr/bin/env bash
# Store the OpenAI runtime API key for tunnel-client in a private file.
# The key is read without echo and never placed on a command line or in shell history.
set -euo pipefail
source "$(dirname "$0")/common.sh"
umask 077
mkdir -p "$CONFIG_DIR"
read -rsp "Paste the tunnel runtime API key (input hidden), then press Enter: " key; echo
[[ "$key" == sk-* ]] || { echo "That does not look like an OpenAI API key (expected sk-...). Nothing saved." >&2; exit 1; }
printf '%s\n' "$key" > "$KEY_FILE"
chmod 600 "$KEY_FILE"
echo "Saved to $KEY_FILE (mode 600)."
