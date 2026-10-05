#!/usr/bin/env bash
# Download the two third-party binaries MolTalk uses, verified against their published SHA-256 checksums:
#   OPSIN 2.9.0 (IUPAC name parser, needs Java)  -> vendor/
#   tunnel-client 0.0.15 (OpenAI Secure MCP Tunnel, Linux x86_64) -> tools/tunnel-client/
set -euo pipefail
cd "$(dirname "$0")/.."

OPSIN=opsin-cli-2.9.0-jar-with-dependencies.jar
OPSIN_SHA=c2e29326c281f87b59a05d934d8589adac6e9d17b95b984931b3e739111b360f
mkdir -p vendor
if [[ ! -f vendor/$OPSIN ]]; then
  curl -fsSL -o vendor/$OPSIN.part "https://github.com/dan2097/opsin/releases/download/2.9.0/$OPSIN"
  echo "$OPSIN_SHA  vendor/$OPSIN.part" | sha256sum -c - && mv vendor/$OPSIN.part vendor/$OPSIN
fi

TC=tunnel-client-v0.0.15-linux-amd64.zip
mkdir -p tools/tunnel-client
if [[ ! -x tools/tunnel-client/tunnel-client ]]; then
  base=https://github.com/openai/tunnel-client/releases/download/v0.0.15
  curl -fsSL -o tools/tunnel-client/$TC "$base/$TC"
  curl -fsSL -o tools/tunnel-client/SHA256SUMS.txt "$base/SHA256SUMS.txt"
  (cd tools/tunnel-client && grep " $TC\$" SHA256SUMS.txt | sha256sum -c - && unzip -oq $TC && rm $TC)
fi
echo "OPSIN: vendor/$OPSIN"
echo "tunnel-client: tools/tunnel-client/tunnel-client"
