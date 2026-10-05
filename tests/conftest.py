import os

# PubChem name lookups only in the network tests (MOLTALK_NETWORK_TESTS=1); everything else runs offline.
if os.getenv("MOLTALK_NETWORK_TESTS") != "1":
    os.environ.setdefault("MOLTALK_OFFLINE", "1")
