# MolTalk privacy policy

*Effective 7 October 2026*

MolTalk is an open-source chemistry tool for AI assistants, provided by Dominic Fico. This policy describes the
public MolTalk service at `https://moltalk-411294000488.us-central1.run.app/mcp`, used through ChatGPT, Claude or
another MCP client. The source code is public at <https://github.com/DomFico/moltalk>, so every statement here can be
checked against it.

## What MolTalk receives

When your assistant uses MolTalk, it sends MolTalk a tool request: typically a compound name, a chemical structure
(SMILES), and display options. MolTalk needs these to compute and draw the structure. MolTalk does not receive your
full conversation, name, email address, account profile, or files unless their contents are explicitly included in a
tool request. The MCP host may also provide limited technical metadata, such as an anonymous per-user identifier,
protocol version, user-agent string, and network address used for security, rate limiting, and operation of the
service.

MolTalk has no accounts, no sign-in, no cookies, no advertising and no analytics or tracking.

## What MolTalk keeps

- **No raw chemistry requests are stored.** Structures and names are processed in memory to answer the request and
  are not written to any database or file.
- **Request logs.** For each request the server logs:
  - the operation name (for example, which tool was called);
  - a short digest of the request arguments. Raw tool arguments are not logged. MolTalk records only the first 8
    hexadecimal characters of a SHA-256 digest of the arguments, used to distinguish repeated calls;
  - the MCP protocol version;
  - the first 40 characters of the client's user-agent string;
  - the response status and how long the tool took.

  For a rejected request (an HTTP 4xx error) the error message is logged too. The hosting platform, Google Cloud Run,
  also keeps standard request logs (time, URL path, status, client IP address and user agent). Logs are kept in
  Google Cloud Logging under its default retention (30 days) and are used only to operate and debug the service.
- **Short-lived memory.**
  - To prevent abuse, the server counts requests per user for one minute, keyed by the anonymous user identifier
    ChatGPT sends (`openai/subject`) or, failing that, the client IP address.
  - To avoid showing the same drawing twice when a host repeats a request, it remembers a SHA-256 digest of each
    chat session's last drawing request for 20 seconds.
  - Results of name lookups may be cached in memory while the server instance runs.

  None of this is written to disk, and it disappears when the server instance stops.

## Who else receives data

- **Google Cloud** hosts the service (Cloud Run, in the United States) and stores the logs described above.
- **PubChem** (U.S. National Center for Biotechnology Information) is contacted in two cases:
  - **Compound names:** when a name is not in MolTalk's built-in library, cannot be read by its offline name parser,
    and the assistant explicitly allows a network lookup, the name is sent to PubChem.
  - **IUPAC names:** to find the IUPAC name of a structure that is not in the built-in library, its InChIKey (a hashed
    identifier of the structure) is sent to PubChem.

  PubChem's own policies apply: <https://www.ncbi.nlm.nih.gov/home/about/policies/>.

MolTalk does not sell or rent personal data and does not intentionally disclose request contents to third parties
except as described above. Your AI-assistant provider (for example, OpenAI or Anthropic) separately processes your
conversation and tool invocation under its own privacy policy. The drawing shown in your chat (the MolTalk viewer)
runs inside your AI assistant's app and makes no network requests of its own.

## Children

MolTalk is not directed to children under 13 and does not knowingly seek personal information from children.
Operational logs may contain the limited technical information described above. Users should not include personal
information in chemistry requests.

## Your choices

Do not put personal or confidential information into chemistry requests. If you run MolTalk yourself from the source
code, set `MOLTALK_OFFLINE=1` to turn off all PubChem lookups.

## Changes and contact

Changes to this policy are published at this address, with a new effective date, and are visible in the repository's
history.

**Privacy requests** (questions about your data, or a request to delete it): email <fico.dominic@gmail.com>. Do not use GitHub Issues
for privacy matters; they are public. For ordinary support, see [Support](SUPPORT.md).
