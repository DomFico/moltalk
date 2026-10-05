#!/usr/bin/env python3
"""Build MolTalk's bundled compound library (moltalk/data/compounds.json.gz).

Sources:
  * Wikidata (CC0): chemical compounds that have a PubChem CID and an English Wikipedia article - a practical
    definition of "commonly discussed". Gives the English label, aliases and how many Wikipedias cover it.
  * scripts/library_seed_names.txt: common teaching names (sugars, cofactors, reagents...) that have no English
    Wikipedia article of their own; each is resolved to a CID with PubChem's name service.
  * PubChem (NCBI): isomeric SMILES, standard InChIKey, title, IUPAC name (OpenEye Lexichem) and leading synonyms.
Verification (an entry is dropped if a check fails):
  * RDKit parses PubChem's SMILES and recomputes the same standard InChIKey PubChem reports.
  * At most MAX_HEAVY heavy atoms (drops proteins and polymers).
IUPAC metadata:
  * OPSIN rebuilds each IUPAC name; the name is marked verified only if that structure equals the compound's.
  * Parent-structure locants are precomputed with the same pipeline drawings use (moltalk.naming).
Names are stored in tiers, most trusted first: "seed" (the seed list), "label" (the compound's own Wikidata label
and PubChem title), "alias" (Wikidata aliases) and "names" (all of these plus PubChem's synonym tail). Lower tiers
contain errors - Wikidata lists "ozone" as an alias of phencyclidine and "LSD" for lysergic acid - so the resolver
uses the most trusted tier that knows the name; a name shared by different structures within that tier is
reported as ambiguous instead of being resolved to one of them.

Run from the project root with the venv:  .venv/bin/python scripts/build_library.py
Needs network (Wikidata, PubChem), Java and vendor/ OPSIN (scripts/fetch-tools.sh). Takes roughly 10-20 minutes.
"""
import gzip
import json
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from rdkit import Chem, RDLogger  # noqa: E402
from moltalk import naming  # noqa: E402

RDLogger.DisableLog("rdApp.*")
OUT = ROOT / "moltalk" / "data" / "compounds.json.gz"
UA = "MolTalk-library-builder/1.0 (https://github.com/DomFico/moltalk)"
MAX_HEAVY = 150
SYNONYMS_PER_COMPOUND = 8
SEEDS = ROOT / "scripts" / "library_seed_names.txt"
REGISTRY = re.compile(r"^(\d{2,7}-\d{2}-\d|\d{3}-\d{3}-\d|RefChem:.*|[A-Z]{2,}[-_ ]?\d+|UNII-.*|CHEBI:\d+|CHEMBL\d+|DTXSID\d+|EINECS.*|NSC[- ]?\d+|"
                      r"SCHEMBL\d+|ZINC\d+|AKOS\d+|MFCD\d+|HSDB.*|BRN .*|CCRIS .*|.*[0-9A-Z]{8,}.*)$")


def log(*args):
    print(time.strftime("%H:%M:%S"), *args, flush=True)


def http_json(url, data=None, tries=6):
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": UA, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=180) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if attempt == tries - 1:
                raise
        except (urllib.error.URLError, TimeoutError, ValueError):  # ValueError: truncated JSON (server-side timeout)
            if attempt == tries - 1:
                raise
        time.sleep(2 ** attempt)


def sparql(query):
    url = "https://query.wikidata.org/sparql?" + urllib.parse.urlencode({"query": query, "format": "json"})
    return http_json(url)["results"]["bindings"]


def wikidata():
    base = "?item wdt:P662 ?cid . ?article schema:about ?item ; schema:isPartOf <https://en.wikipedia.org/> ."
    items = {}
    for row in sparql(f"SELECT ?item ?cid ?label ?links WHERE {{ {base} ?item wikibase:sitelinks ?links . "
                      f"OPTIONAL {{ ?item rdfs:label ?label FILTER(LANG(?label) = 'en') }} }}"):
        qid = row["item"]["value"].rsplit("/", 1)[1]
        entry = items.setdefault(qid, {"cids": set(), "label": None, "aliases": set(), "links": 0})
        if row["cid"]["value"].isdigit():
            entry["cids"].add(int(row["cid"]["value"]))
        entry["label"] = entry["label"] or (row.get("label") or {}).get("value")
        entry["links"] = max(entry["links"], int(row["links"]["value"]))
    for row in sparql(f"SELECT ?item ?alias WHERE {{ {base} ?item skos:altLabel ?alias FILTER(LANG(?alias) = 'en') }}"):
        qid = row["item"]["value"].rsplit("/", 1)[1]
        if qid in items:
            items[qid]["aliases"].add(row["alias"]["value"])
    return items


class PubChem:
    def __init__(self):
        self.last = 0.0

    def post(self, path, cids):
        wait = self.last + 0.25 - time.monotonic()  # stay under PubChem's 5 requests per second
        if wait > 0:
            time.sleep(wait)
        self.last = time.monotonic()
        data = urllib.parse.urlencode({"cid": ",".join(map(str, cids))}).encode()
        return http_json("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + path, data=data)

    def cid_for_name(self, name):
        wait = self.last + 0.25 - time.monotonic()
        if wait > 0:
            time.sleep(wait)
        self.last = time.monotonic()
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + urllib.parse.quote(name, safe="") + "/cids/JSON"
        for attempt in range(6):  # "server busy" is not "not found": wait and ask again
            try:
                got = http_json(url, tries=1)
            except urllib.error.HTTPError as e:
                if e.code != 503:
                    raise
                time.sleep(5 * (attempt + 1))
                continue
            cids = (got or {}).get("IdentifierList", {}).get("CID", [])
            return cids[0] if cids and cids[0] else None
        raise RuntimeError(f"PubChem stayed busy for seed name {name!r}")

    def properties(self, cids):
        out = {}
        for i in range(0, len(cids), 200):
            got = self.post("property/SMILES,InChIKey,Title,IUPACName,HeavyAtomCount/JSON", cids[i:i + 200])
            for row in (got or {}).get("PropertyTable", {}).get("Properties", []):
                out[row["CID"]] = row
            if i % 4000 == 0:
                log(f"  properties {i + 200}/{len(cids)}")
        return out

    def synonyms(self, cids):
        out = {}
        for i in range(0, len(cids), 100):
            got = self.post("synonyms/JSON", cids[i:i + 100])
            for row in (got or {}).get("InformationList", {}).get("Information", []):
                names = [s for s in row.get("Synonym", []) if not REGISTRY.match(s) and len(s) <= 80]
                out[row["CID"]] = names[:SYNONYMS_PER_COMPOUND]
            if i % 4000 == 0:
                log(f"  synonyms {i + 100}/{len(cids)}")
        return out


def opsin_batch(names):
    """{name: single-molecule CML string or None}, one JVM per chunk."""
    jar = naming.OPSIN_JAR
    out = {}
    for i in range(0, len(names), 2000):
        chunk = names[i:i + 2000]
        proc = subprocess.run(["java", "-Xmx1g", "-jar", str(jar), "-ocml"], input="\n".join(chunk) + "\n",
                              capture_output=True, text=True, timeout=1800)
        root = ET.fromstring(proc.stdout[proc.stdout.index("<?xml"):]) if "<?xml" in proc.stdout else None
        molecules = list(root.iter(naming.CML + "molecule")) if root is not None else []
        if len(molecules) != len(chunk):
            raise RuntimeError(f"OPSIN returned {len(molecules)} results for {len(chunk)} names")
        for name, molecule in zip(chunk, molecules):
            out[name] = ET.tostring(molecule, encoding="unicode") if molecule.find(naming.CML + "atomArray") is not None else None
        log(f"  OPSIN {i + len(chunk)}/{len(names)}")
    return out


def normalise(name: str) -> str:
    name = name.replace("‐", "-").replace("‑", "-").replace("–", "-").replace("−", "-").replace("’", "'")
    return re.sub(r"\s+", " ", name.strip().lower())


def main():
    log("Wikidata: compounds with a PubChem CID and an English Wikipedia article")
    items = wikidata()
    log(f"  {len(items)} items")
    # One CID per item: the lowest (PubChem's parent record is usually the first registered).
    cid_items = {}
    for qid, entry in items.items():
        if entry["cids"]:
            cid_items.setdefault(min(entry["cids"]), []).append(qid)
    pubchem = PubChem()
    seeds = defaultdict(set)
    seed_names = [l.strip() for l in SEEDS.read_text().splitlines() if l.strip() and not l.startswith("#")]
    log(f"PubChem CIDs for {len(seed_names)} seed names")
    missing = seed_names
    for attempt in range(3):  # PubChem's name service sometimes answers 404 for names it knows: ask again later
        if attempt:
            log(f"  {len(missing)} seed names not found; asking again in 60 s")
            time.sleep(60)
        retry = []
        for line in missing:
            name, _, pinned = (part.strip() for part in line.partition("|"))
            cid = int(pinned) if pinned else pubchem.cid_for_name(name)
            if cid:
                seeds[cid].add(name)
            else:
                retry.append(line)
        missing = retry
        if not missing:
            break
    if missing:
        raise RuntimeError(f"seed names PubChem does not know (fix {SEEDS.name}): {missing}")
    cids = sorted(set(cid_items) | set(seeds))
    log(f"PubChem properties for {len(cids)} CIDs")
    props = pubchem.properties(cids)
    log(f"PubChem synonyms for {len(props)} CIDs")
    synonyms = pubchem.synonyms(sorted(props))

    log("RDKit verification")
    compounds, dropped = [], defaultdict(int)
    for cid in sorted(props):
        row = props[cid]
        smiles, key = row.get("SMILES"), row.get("InChIKey")
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        if mol is None:
            dropped["unparsable"] += 1
            continue
        if mol.GetNumAtoms() > MAX_HEAVY:
            dropped["too large"] += 1
            continue
        if Chem.MolToInchiKey(mol) != key:
            dropped["InChIKey mismatch"] += 1
            continue
        seed, label, alias = set(seeds.get(cid, ())), set(), set()
        links = 0
        for qid in cid_items.get(cid, ()):
            item = items[qid]
            links = max(links, item["links"])
            label |= {item["label"]} if item["label"] else set()
            alias |= set(item["aliases"])
        if row.get("Title"):
            label.add(row["Title"])
        names = seed | label | alias | set(synonyms.get(cid, []))
        compounds.append({"cid": cid, "title": row.get("Title") or (sorted(names)[0] if names else str(cid)),
                          "smiles": smiles, "inchikey": key, "iupac": row.get("IUPACName"), "links": links,
                          "seed": sorted(seed), "label": sorted(label),
                          "alias": sorted({n for n in alias if 0 < len(n) <= 120}),
                          "names": sorted({n for n in names if 0 < len(n) <= 120})})
    log(f"  kept {len(compounds)}; dropped {dict(dropped)}")

    log("OPSIN verification of IUPAC names and parent locants")
    parsed = opsin_batch(sorted({c["iupac"] for c in compounds if c["iupac"]}))
    stats = defaultdict(int)
    for n, c in enumerate(compounds):
        cml = parsed.get(c["iupac"]) if c["iupac"] else None
        c["iupac_verified"] = bool(cml) and naming.matches_structure(c["smiles"], cml)
        if c["iupac_verified"]:
            try:
                locants, _ = naming.locants_from_cml(c["smiles"], c["iupac"], cml)
            except Exception:
                locants = None
            if locants:
                c["locants"] = locants
                stats["locants"] += 1
            stats["verified"] += 1
        if n % 2000 == 0:
            log(f"  {n}/{len(compounds)}")
    log(f"  IUPAC names verified by OPSIN: {stats['verified']}; with parent locants: {stats['locants']}")

    tiers = [defaultdict(set) for _ in range(4)]
    for i, c in enumerate(compounds):
        iupac = [c["iupac"]] if c["iupac"] else []
        for tier, names in zip(tiers, (c["seed"], c["label"] + iupac, c["alias"], c["names"] + iupac)):
            for name in names:
                tier[normalise(name)].add(i)
    index = tiers[-1]
    ambiguous = sum(1 for n in index
                    if len({compounds[i]["inchikey"] for i in next(t[n] for t in tiers if n in t)}) > 1)
    log(f"  {len(index)} distinct names ({ambiguous} still shared by different structures: reported as ambiguous)")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    payload = {"version": 1, "built": date.today().isoformat(),
               "sources": {"names": "Wikidata (CC0): compounds with a PubChem CID and an English Wikipedia article; "
                                    "plus scripts/library_seed_names.txt resolved with PubChem",
                           "structures": "PubChem (NCBI): isomeric SMILES, InChIKey, title, IUPAC name, synonyms",
                           "verification": "RDKit InChIKey re-derivation; OPSIN name-to-structure for IUPAC names"},
               "compounds": compounds}
    with gzip.open(OUT, "wt", encoding="utf-8") as f:
        json.dump(payload, f, separators=(",", ":"))
    log(f"Wrote {OUT} ({OUT.stat().st_size / 1e6:.1f} MB, {len(compounds)} compounds)")


if __name__ == "__main__":
    main()
