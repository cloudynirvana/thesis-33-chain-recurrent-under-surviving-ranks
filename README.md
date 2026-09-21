# Chain-recurrent labels recoverable from multi-channel ranks that survive stiff–sloppy reduction

**Thesis #33.** Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Depends on:** Thesis #26 (chain-recurrent versus hybrid occult) and Thesis #24 (reduction-preserving multi-channel identifiability).

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

Which chain-recurrent versus hybrid-occult labels remain recoverable when observation is restricted to the multi-channel Fisher ranks that Thesis #24 proved survive a stiff–sloppy reduction?

On a joint hybrid–smooth toy they split by label. The proliferative field has an unstable focus and a periodic orbit of period 19.874. Two occult guards are read on that field. In 240 stratified trials every trajectory hits a guard, so the occult class is a two-way split and not a three-way that includes a no-hit class. A multi-channel map with a stiff modifier has practical Fisher rank 3 of 5; after the documented reduction that keeps κ = h a / b the reduced practical rank is 3 of 3, and those three full directions survive. A linear reader sits at chance on the cycle label under the full feature and under the surviving-rank projection (hold-out accuracies 0.404 and 0.392 against chance 0.5). The occult-guard label remains recoverable from the surviving ranks (hold-out accuracy 1.000). The orbital check is positive on 54.0 percent of the orbit; that is a collocation defect, not a certified Conley set.

No number is taken from either parent deposit's results file. This deposit does not claim a staging rule or a cure.

This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md).

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `sim/surviving_ranks_labels.py` | Joint toy: guards, Fisher ranks, linear readers (seed 20260921) |
| `sim/results.json` | Numbers cited in Chapter Four |
| `sim/figures/` | Label recovery, cycle cloud, guard hits, Fisher spectra |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/surviving_ranks_labels.py
python3 build_pdf.py
```

NumPy, SciPy and Matplotlib are required for the toy. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json` and `sim/figures/`.

## Cite

Ogbonna KE. Chain-recurrent labels recoverable from multi-channel ranks that survive stiff–sloppy reduction [Internet]. Thesis #33 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-33-chain-recurrent-under-surviving-ranks

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.
