"""Public HTTP gateway: health check, domain verification, per-user rate limits."""
import importlib
from starlette.testclient import TestClient

HEADERS = {"Accept": "application/json, text/event-stream"}


def _client(monkeypatch, **env):
    for k, v in env.items():
        monkeypatch.setenv(k, v)
    from moltalk import public, server
    importlib.reload(public)
    server.mcp._session_manager = None
    return TestClient(public.PublicGateway(server.mcp.streamable_http_app()), base_url="http://localhost:8000")


def _list(client, subject):
    return client.post("/mcp", headers=HEADERS, json={"jsonrpc": "2.0", "id": 1, "method": "tools/list",
                                                      "params": {"_meta": {"openai/subject": subject}}}).status_code


def test_health_and_domain_verification(monkeypatch):
    with _client(monkeypatch, MOLTALK_OPENAI_CHALLENGE="token-abc") as client:
        assert client.get("/health").text == "ok"
        assert client.get("/.well-known/openai-apps-challenge").text == "token-abc"


def test_rate_limit_is_per_chatgpt_user(monkeypatch):
    with _client(monkeypatch, MOLTALK_RATE_PER_USER="3") as client:
        assert [_list(client, "alice") for _ in range(4)] == [200, 200, 200, 429]
        assert _list(client, "bob") == 200  # other users share OpenAI's IPs but are limited separately
