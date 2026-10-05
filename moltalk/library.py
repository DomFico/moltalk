"""MolTalk's bundled compound library (built by scripts/build_library.py; see that script for sources).

Thousands of commonly discussed compounds with PubChem CID, isomeric SMILES, standard InChIKey, title, synonyms,
PubChem's IUPAC name (marked verified when OPSIN rebuilds exactly that structure) and precomputed parent locants.
Lookups are local: no network, nothing sent anywhere.
"""
import gzip
import json
import re
from collections import defaultdict
from functools import lru_cache
from importlib.resources import files

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

DATA = files("moltalk").joinpath("data/compounds.json.gz")
# A stereo qualifier in front of a name: L-, D-, DL-, rac-, meso-, cis-, trans-, (+)-, (-)-, (±)-, (R)-, (2S,3R)-, (E)-...
_QUALIFIER = r"(?:[dl]|dl|d,l|rac|meso|cis|trans|[+\-±]|\([+\-±]+\)|\([0-9a-z,'\s]*[rsez]\)|[dl]-\([+\-]\)|\([+\-]\)-[dl])"
_QUALIFIED = re.compile(r"^" + _QUALIFIER + r"[- ]")
# Qualifiers that make a bare name a stereo family (amino acids, sugars): D/L (Fischer) and racemic. (+)/(-), R/S and
# cis/trans are not included: "morphine" names one enantiomer, and ephedrine/pseudoephedrine are different compounds.
_FAMILY = r"(?:[dl]|dl|d,l|rac|±|\(±\)|\(\+-\)|\(\+/-\)|[dl]-\([+\-]\)|\([+\-]\)-[dl])"


def normalise(name: str) -> str:
    """Case-, whitespace- and dash-insensitive form used for name lookups."""
    name = name.replace("‐", "-").replace("‑", "-").replace("–", "-").replace("−", "-").replace("’", "'")
    return re.sub(r"\s+", " ", name.strip().lower())


@lru_cache(maxsize=1)
def _load():
    if not DATA.is_file():
        return {"compounds": [], "tiers": [{}], "by_key": {}, "by_skeleton": {}, "built": None}
    payload = json.loads(gzip.decompress(DATA.read_bytes()))
    compounds = payload["compounds"]
    tiers, by_key, by_skeleton = [defaultdict(list) for _ in range(4)], {}, defaultdict(list)
    for i, c in enumerate(compounds):
        iupac = [c["iupac"]] if c.get("iupac") else []
        for tier, names in zip(tiers, (c.get("seed", []), c.get("label", []) + iupac, c.get("alias", []),
                                       c["names"] + iupac)):
            for name in names:
                tier[normalise(name)].append(i)
        by_skeleton[c["inchikey"][:14]].append(i)
        best = by_key.get(c["inchikey"])
        if best is None or c["cid"] < compounds[best]["cid"]:  # one record per structure: PubChem's lowest CID
            by_key[c["inchikey"]] = i
    return {"compounds": compounds, "tiers": [dict(t) for t in tiers], "by_key": by_key,
            "by_skeleton": dict(by_skeleton), "built": payload.get("built")}


def size() -> int:
    return len(_load()["compounds"])


def built() -> str | None:
    return _load()["built"]


def _has_stereo(smiles: str) -> bool:
    return any(ch in smiles for ch in "@/\\")


def stereo_alternatives(name: str, compound: dict):
    """Other stereoisomers in the library that the same unqualified name could mean, or None if the name is fine.
    "alanine" -> L-alanine is a convention, not something the name says. A bare name is treated as a stereo family
    when the record's own title/label carries a D/L or racemic qualifier (L-Alanine, D-Fructose), or when another
    stereoisomer's own title/label is that qualifier plus the name (D-Alanine, DL-alanine). Synonyms do not count:
    PubChem lists "d-morphine" for (+)-morphine, but "morphine" names the natural enantiomer.
    None when the name is itself stereo-qualified, the structure has no stereo, or the entry is curated."""
    key = normalise(name)
    if compound.get("curated") or _QUALIFIED.match(key) or not _has_stereo(compound["smiles"]):
        return None
    data = _load()
    pattern = re.compile(r"^" + _FAMILY + r"[- ]" + re.escape(key) + r"$")
    own = lambda c: [normalise(n) for n in c.get("seed", []) + c.get("label", []) + [c["title"]]]
    found = {}
    for i in data["by_skeleton"].get(compound["inchikey"][:14], []):
        other = data["compounds"][i]
        if other["inchikey"] != compound["inchikey"] and other["inchikey"] not in found \
                and any(pattern.match(n) for n in own(other)):
            found[other["inchikey"]] = other
    if found or any(pattern.match(n) for n in own(compound)):
        return sorted(found.values(), key=lambda c: c["cid"])
    return None


def find_name(name: str):
    """("hit", compound) | ("ambiguous", [compounds]) | ("stereo", [compounds]) | ("miss", None).
    "stereo": the name matches one record, but it does not say which stereoisomer, and the library holds others
    under that name with a stereo qualifier (see stereo_alternatives); the hit comes first.
    Name tiers, most trusted first: seed list; the compound's own Wikidata label, PubChem title or IUPAC name;
    Wikidata aliases; PubChem synonyms. The first tier that knows the name decides. A name shared by different
    structures within that tier is reported as ambiguous rather than resolved to one of them."""
    data = _load()
    key = normalise(name)
    hits = next((t[key] for t in data["tiers"] if key in t), [])
    unique = {}
    for i in hits:
        c = data["compounds"][i]
        if c["inchikey"] not in unique or c["cid"] < unique[c["inchikey"]]["cid"]:
            unique[c["inchikey"]] = c
    if not unique:
        return "miss", None
    if len(unique) == 1:
        hit = next(iter(unique.values()))
        others = stereo_alternatives(name, hit)
        return ("stereo", [hit, *others]) if others is not None else ("hit", hit)
    return "ambiguous", sorted(unique.values(), key=lambda c: -c.get("links", 0))


def by_inchikey(key: str):
    data = _load()
    i = data["by_key"].get(key)
    return data["compounds"][i] if i is not None else None


def locants_for(smiles: str, compound: dict):
    """Map the library's precomputed parent locants (indexed on its own SMILES) onto the user's molecule.
    Same rules as live numbering: any symmetry-equivalent numbering is valid, but only one under which every
    R/S descriptor in the IUPAC name matches the molecule's own CIP label. ({index: locant}, None) or (None, reason)."""
    from .naming import _map_locants, _stereo_descriptors
    if not compound.get("locants"):
        return None, ("the IUPAC name could not be verified with OPSIN" if not compound.get("iupac_verified")
                      else "the parent chain or ring of the name could not be identified unambiguously")
    mol = Chem.MolFromSmiles(smiles)
    ref = Chem.MolFromSmiles(compound["smiles"])
    rdCIPLabeler.AssignCIPLabels(mol)
    parent = {int(k): v for k, v in compound["locants"].items()}
    found = _map_locants(mol, ref, parent, _stereo_descriptors(compound["iupac"]))
    if not found:
        return None, "no numbering of the parent matches the name's stereodescriptors"
    return {str(i): l for i, l in sorted(found.items())}, None
