# qPCR analysis in R — tutorial

Data and code to follow along with a walkthrough of basic qPCR analysis in R.

**Full walkthrough:** <https://www.patrickfhorve.com/bioinformatics/qpcr>

## What it covers

Starting from raw instrument output, the walkthrough works through:

- reading an Applied Biosystems export into R and handling `Undetermined` wells
- averaging technical replicates by sample
- separating standards from unknowns
- building a standard curve by linear regression of Ct against log10 gene copy number
- plotting the curve with `ggplot2` and reading the fit back off it
- converting sample Ct values into gene copy numbers

## Contents

| File | What it is |
|---|---|
| `qPCR_Tutorial.R` | The analysis script the walkthrough builds up |
| `PMATest_data.csv` | Example dataset — a SYBR Green run on a 96-well block, exported from an Applied Biosystems instrument |

## Requirements

R, with `dplyr` and `ggplot2`.

## Use

Clone the repository, point the `read.csv()` call at `PMATest_data.csv`, and work through
`qPCR_Tutorial.R` alongside the walkthrough linked above.

## License

MIT — see [LICENSE](LICENSE).
