# Data sources

| Source | Used for | Licence |
|---|---|---|
| [RDKit](https://www.rdkit.org/) | Parsing, validation, CIP labels, 2D layout and drawing, 3D embedding and force fields, file formats | BSD 3-Clause |
| [OPSIN](https://github.com/dan2097/opsin) 2.9.0 | Reading systematic IUPAC names; checking library and PubChem names; locants | MIT (downloaded, not in the repository) |
| [PubChem](https://pubchem.ncbi.nlm.nih.gov/) | Library structures and names; optional live fallback | NCBI [policies](https://www.ncbi.nlm.nih.gov/home/about/policies/) |
| [Wikidata](https://www.wikidata.org/) | Which compounds are in the library, and their labels and aliases | CC0 |
| [MCP Python SDK](https://github.com/modelcontextprotocol/python-sdk) | The MCP server | MIT |

## The bundled library

`moltalk/data/compounds.json.gz` (about 3 MB) holds 22,038 compounds: every compound on Wikidata that has a PubChem
CID and an English Wikipedia article, plus about 460 common teaching names from `scripts/library_seed_names.txt`
(sugars, cofactors, terpenes, common drugs, reagents). It covers drugs, natural products, metabolites, amino acids,
solvents and reagents, with about 200,000 names.

Each entry has its PubChem CID, PubChem's isomeric SMILES and title, synonyms, PubChem's IUPAC name, whether OPSIN
rebuilds exactly that structure from the name (about 89 %), and precomputed parent locants (about 60 %).

**Build rules** (`scripts/build_library.py`): RDKit must reproduce PubChem's InChIKey from its SMILES, or the entry is
dropped; compounds over 150 heavy atoms are dropped; any seed name PubChem cannot find stops the build.
`scripts/library_curated.json` holds editorial overrides (heme vs hemin, NAD⁺/NADH and other cofactors), each with a
corrected title and a note returned with the structure.

## Live lookups

- **Name → structure** reaches PubChem only when the name is not in the library, OPSIN cannot read it, and the model
  passes `allow_network=true`; the name is then sent to PubChem.
- **Structure → name**: a structure not in the library has its InChIKey looked up on PubChem (cached). Set
  `MOLTALK_OFFLINE=1` to disable all PubChem access; library names and numbering keep working.
