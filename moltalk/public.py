"""Public HTTP deployment (Cloud Run): rate limits, health check, OpenAI domain verification.

Requests from ChatGPT come from OpenAI's servers, so many users share a few IP addresses. Limits are therefore
keyed by ChatGPT's anonymous per-user id (params._meta["openai/subject"]) when present, and by client IP otherwise,
with a global cap on top so a traffic spike cannot run up the hosting bill.
"""
import json
import os
import time

PER_USER_PER_MINUTE = int(os.getenv("MOLTALK_RATE_PER_USER", "60"))
GLOBAL_PER_MINUTE = int(os.getenv("MOLTALK_RATE_GLOBAL", "600"))
CHALLENGE = os.getenv("MOLTALK_OPENAI_CHALLENGE", "").strip()


class _Window:
    """Requests per key in the last 60 seconds (sliding window)."""

    def __init__(self, limit: int):
        self.limit, self.hits = limit, {}

    def allow(self, key: str, now: float) -> bool:
        recent = [t for t in self.hits.get(key, ()) if now - t < 60]
        if len(recent) >= self.limit:
            self.hits[key] = recent
            return False
        recent.append(now)
        self.hits[key] = recent
        if len(self.hits) > 50_000:  # forget idle keys
            self.hits = {k: v for k, v in self.hits.items() if v and now - v[-1] < 60}
        return True


def _client_ip(scope) -> str:
    for name, value in scope.get("headers", []):
        if name == b"x-forwarded-for":  # set by Cloud Run's front end; the first address is the client
            return value.decode(errors="replace").split(",")[0].strip()
    client = scope.get("client")
    return client[0] if client else "unknown"


def _subject(body: bytes) -> str | None:
    try:
        message = json.loads(body)
    except (ValueError, UnicodeDecodeError):
        return None
    messages = message if isinstance(message, list) else [message]
    for m in messages:
        meta = ((m or {}).get("params") or {}).get("_meta") or {}
        if isinstance(meta, dict) and meta.get("openai/subject"):
            return "user:" + str(meta["openai/subject"])[:200]
    return None


async def _send_plain(send, status: int, text: str, content_type: bytes = b"text/plain; charset=utf-8"):
    body = text.encode()
    await send({"type": "http.response.start", "status": status,
                "headers": [(b"content-type", content_type), (b"content-length", str(len(body)).encode())]})
    await send({"type": "http.response.body", "body": body})


class PublicGateway:
    """ASGI wrapper around the MCP app."""

    def __init__(self, app):
        self.app = app
        self.per_user, self.overall = _Window(PER_USER_PER_MINUTE), _Window(GLOBAL_PER_MINUTE)

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http":
            return await self.app(scope, receive, send)
        path = scope.get("path", "")
        if path == "/healthz":
            return await _send_plain(send, 200, "ok")
        if path == "/.well-known/openai-apps-challenge":
            return await _send_plain(send, 200, CHALLENGE) if CHALLENGE else await _send_plain(send, 404, "not configured")
        if scope.get("method") != "POST":
            return await self.app(scope, receive, send)

        # Read the body once (it is small: the MCP server caps request size), then replay it to the MCP app.
        chunks, more = [], True
        while more:
            message = await receive()
            chunks.append(message.get("body", b""))
            more = message.get("more_body", False)
        body = b"".join(chunks)
        now = time.monotonic()
        key = _subject(body) or "ip:" + _client_ip(scope)
        if not self.overall.allow("all", now) or not self.per_user.allow(key, now):
            error = {"jsonrpc": "2.0", "id": None,
                     "error": {"code": -32000, "message": "MolTalk is busy or you are sending requests too quickly. "
                                                          "Please wait a minute and try again."}}
            return await _send_plain(send, 429, json.dumps(error), b"application/json")

        replayed = False

        async def replay():
            nonlocal replayed
            if not replayed:
                replayed = True
                return {"type": "http.request", "body": body, "more_body": False}
            return await receive()
        return await self.app(scope, replay, send)
