"""MolTalk: chemistry tools for AI assistants (RDKit, OPSIN and a curated compound library over MCP)."""
from importlib.metadata import PackageNotFoundError, version as _version

try:
    __version__ = _version("moltalk")
except PackageNotFoundError:  # running from a source tree without installation
    __version__ = "0+unknown"
