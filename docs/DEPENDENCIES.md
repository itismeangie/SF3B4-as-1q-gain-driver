# Dependencies

This repository contains computational scripts only. It does not vendor package libraries, raw data, generated outputs, or local caches.

## Runtime

The manuscript software table should be treated as the authoritative version record at final submission. The current public-code package expects:

- R 4.4 or newer. The manuscript environment used R 4.5.1.
- Python 3.10 or newer. The manuscript environment used Python 3.13.5.
- A POSIX shell for `run.sh` wrappers.
- `curl` for selected public downloads.

## Python Packages

The Python scripts use:

- `numpy`
- `openpyxl`
- `pandas`

Install with:

```bash
python3 -m pip install -r requirements.txt
```

## R Packages

The R scripts use CRAN packages:

- `broom`
- `data.table`
- `dplyr`
- `readr`
- `readxl`
- `survival`
- `tibble`
- `tidyr`
- `tidyverse`

They also use Bioconductor packages:

- `Biobase`
- `BiocParallel`
- `DESeq2`
- `edgeR`
- `fgsea`
- `limma`

NetBID2 analyses additionally require `NetBID2`, which should be installed according to the NetBID2 project instructions used for the manuscript environment.

The SJAracne network runner also expects `sjaracne` on `PATH` when `04_netbid2_driver_prioritization/run_sjaracne_networks.sh` is used.

## Conda Skeleton

`environment.yml` provides a lightweight starting point for common dependencies. It may still require manual installation of `NetBID2` and any site-specific packages used with controlled local data.

```bash
conda env create -f environment.yml
conda activate sf3b4-1q-gain-driver
```

## Controlled Inputs

Dependencies alone are not sufficient to rerun all topics. Several analyses require controlled or local external inputs listed in `docs/EXTERNAL_INPUTS.md`.
