# MolTalk support

## Getting help

- **Report a problem or ask a question:** open an issue at <https://github.com/DomFico/moltalk/issues>. Include the
  prompt you used, what you expected, and what happened. A screenshot of the drawing helps. Please do not include
  personal or confidential information.
- **Documentation:** start with the [README](../README.md); details are in [Chemistry](CHEMISTRY.md),
  [Scope and limitations](LIMITATIONS.md) and [MCP tools and the viewer](MCP_AND_UI.md).

## Common questions

**The drawing does not appear, only text.** Refresh the MolTalk app or connector in your assistant and start a new
chat; hosts cache the viewer. The assistant still receives the full chemistry data.

**"No structure was assumed" or "does not say which stereoisomer".** MolTalk refuses names it cannot resolve to
exactly one structure, instead of guessing. Give a more specific name (for example "D-glucose"), a SMILES string, or
ask the assistant to list the candidates.

**A name with (R)/(S) or (P)/(M) is refused** (for example "(R)-BINAP"). The databases MolTalk uses record such
compounds without their axial or helical configuration, so it cannot verify the descriptor. Ask for "the
stereoisomers of BINAP" to get both forms with verified descriptors.

**The service says it is busy.** There is a per-user limit of 60 requests per minute; wait a minute and try again.

**I want to run MolTalk myself.** It is open source: see [Run it yourself](../README.md#run-it-yourself) and
[Deployment](DEPLOYMENT.md).

## Policies

[Privacy policy](PRIVACY.md) · [Terms of service](TERMS.md)
