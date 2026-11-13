# mscode

R package collecting helper functions for MS-based proteomics data analysis at different
levels: ID overlap and set comparisons, correlation networks and benchmarking, sibling-peptide
correlation, FASTA/ID parsing, format conversion (mzIdentML, FragPipe, MaxLynx → xiVIEW,
IntAct), and QC plots (TIC profiles, peptide profiles, filter impact, PANTHER results).

## Installation

```r
devtools::install_github("heuselm/mscode")
library(mscode)
```

## Standalone scripts

`scripts/` holds analysis scripts that are not part of the package build:

| Script | Purpose |
|---|---|
| [`processObisDDAresults/`](scripts/processObisDDAresults) | Filter, summarize and visualize TPP/iProphet DDA-MS identification results at a target FDR using the Mayu FDR table (formerly the separate `processObisDDAresults` repo; history preserved) |

## License

MIT — see [LICENSE](LICENSE).

---
<sub>README optimized with Claude Opus 5.5.</sub>
