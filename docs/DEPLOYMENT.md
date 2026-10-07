# Deployment

MolTalk's public instance runs on **Google Cloud Run** at

```text
https://moltalk-411294000488.us-central1.run.app/mcp
```

It is stateless, scales to zero, and needs nobody's computer to be on.

## Cloud Run (production)

From the repository root:

```bash
gcloud run deploy moltalk --source . --region us-central1 --allow-unauthenticated \
  --min-instances 0 --max-instances 2 --cpu 1 --memory 2Gi --concurrency 8 --timeout 60 --cpu-boost \
  --set-env-vars "^|^MOLTALK_ALLOWED_HOSTS=moltalk-411294000488.us-central1.run.app,moltalk-jo5jxgfila-uc.a.run.app,localhost:*,127.0.0.1:*"
```

For later deploys of the same service, `gcloud run deploy moltalk --source . --region us-central1 --quiet` keeps the
existing settings.

- **Image:** the `Dockerfile` installs RDKit's drawing libraries and Java, downloads OPSIN 2.9.0 (SHA-256 verified),
  installs the package and runs as a non-root user. The bundled compound library ships inside the package.
- **Cost:** it scales to zero, so an idle service costs nothing, and ordinary use stays inside Cloud Run's free tier.
  `--max-instances 2` caps the bill.
- **Protection** (`moltalk/public.py`):
  - Host/Origin validation: only the hosts in `MOLTALK_ALLOWED_HOSTS` are accepted.
  - Rate limits: 60 requests per minute per user (ChatGPT's `openai/subject`, or the client IP) and 600 per minute
    in total.
  - A 64 KB request-body limit.
- **Endpoints:** `/mcp` (MCP), `/health` (health check; Cloud Run reserves `/healthz`), `/logo.png` and
  `/favicon.ico`, and `/.well-known/openai-apps-challenge` (serves `MOLTALK_OPENAI_CHALLENGE` for OpenAI's domain
  verification).
- **Logs:** each request's JSON-RPC method and tool name (never the arguments), the client's
  `mcp-protocol-version`, and the reason for any 4xx response. Tool run times appear as `moltalk <tool> took N s`.

### Environment variables

| Variable | Default | Purpose |
|---|---|---|
| `MOLTALK_ALLOWED_HOSTS` / `MOLTALK_ALLOWED_ORIGINS` | localhost only | Accepted Host and Origin headers |
| `MOLTALK_WORKERS` | 2 | RDKit worker processes |
| `MOLTALK_TIMEOUT_S` | 20 (25 in the image) | Per-call timeout; the worker is killed when exceeded |
| `MOLTALK_WORKER_MEMORY_MB` | 2048 | Address-space cap per worker |
| `MOLTALK_MAX_QUEUED` | 16 | Queued requests before "busy" |
| `MOLTALK_MAX_BODY_BYTES` | 65536 | Request-body limit |
| `MOLTALK_OFFLINE` | unset | `1` disables PubChem lookups entirely |
| `MOLTALK_OPSIN_JAR` | `vendor/…` | OPSIN jar path (set in the image) |
| `MOLTALK_OPENAI_CHALLENGE` | unset | OpenAI domain-verification token |

### PubChem and outgoing IPs

The library and OPSIN run inside the container; PubChem is only a fallback, rate-limited to 4 requests per second.
From Cloud Run's shared egress addresses PubChem often answers "server busy". If real users repeatedly need names
that are not in the library, the cheapest fix is to add them to `scripts/library_seed_names.txt` and rebuild the
library. Routing egress through a static IP (VPC connector or Direct VPC egress, Cloud Router, Cloud NAT and a
reserved address, roughly US$4–5 a month) is documented but not provisioned.

## Connecting clients

- **ChatGPT:** the plugin package in `plugin/` points at the Cloud Run URL (see [Plugin package](#plugin-package)).
  In developer mode, a connector can also be added directly with the MCP URL and no authentication.
- **Claude:** add a custom connector with the MCP URL (Settings → Connectors).
- **Other MCP clients:** connect to the URL with the streamable-HTTP transport, or run the server locally over stdio
  (`moltalk`).

After the widget changes, hosts may cache the old template: refresh the plugin or connector and start a new chat.

## Plugin package

`plugin/` is the portable plugin package:

- `plugin.json`: name, descriptions, `logo` and `composerIcon` under `extensions.com.openai.interface`;
- `mcp.json`: the Cloud Run MCP URL;
- `assets/moltalk-logo.png`: a 512×512 PNG, the canonical app icon.

Build the upload with:

```bash
cd plugin && zip -r ../dist/moltalk-plugin.zip .
```

Logo and metadata change only when the package is re-uploaded; redeploying the server does not update them. The
server also advertises its icon in `serverInfo.icons` and serves it at `/logo.png`.

## Docker locally

```bash
docker build -t moltalk .
docker run --rm -p 127.0.0.1:8000:8000 -e PORT=8000 --read-only --tmpfs /tmp --memory 2g --cpus 2 --pids-limit 128 moltalk
```

The server listens on `0.0.0.0` inside the container, but Host/Origin validation still accepts only localhost
unless `MOLTALK_ALLOWED_HOSTS` is set. There is no authentication, so keep a local server on loopback.

## Legacy: OpenAI Secure MCP Tunnel

Before the Cloud Run deployment, MolTalk ran on a personal computer and reached ChatGPT through OpenAI's Secure MCP
Tunnel (outbound-only, no open ports). The scripts are kept for anyone who wants a private, local instance:

```bash
scripts/fetch-tools.sh                    # tunnel-client and OPSIN, checksum-verified
scripts/store-tunnel-key.sh               # stores a Platform runtime key (mode 600); never paste keys into a chat
scripts/configure-tunnel.sh tunnel_…      # writes the tunnel profile and runs tunnel-client doctor
scripts/start.sh                          # systemd user service moltalk-tunnel; --foreground to run in the terminal
scripts/status.sh                         # health; admin UI on 127.0.0.1:8791 only
scripts/stop.sh / scripts/uninstall-service.sh
scripts/local-tunnel-test.sh              # the tunnel path locally, without credentials
```

It needs a ChatGPT plan with developer mode, an OpenAI Platform organization with a tunnel (Tunnels Read + Manage to
create, Read + Use to run), and the computer to stay on while the plugin is used. Create the plugin with
**Connection: Tunnel** and **No authentication**. Logs: `journalctl --user -u moltalk-tunnel -f`.
