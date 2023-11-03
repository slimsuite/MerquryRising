# MerquryRising: Merqury parsing and comparison for genome assembly purging QC

[MerquryRising](https://github.com/slimsuite/MerquryRising) is a standalone RMarkdown file for parsing outputs from [Merqury](https://github.com/marbl/merqury) and generating some additional plots to help assess the requirements and/or consequences of duplicate purging in genome assemblies. Please see the accompanying <Tutorial.html> for additional information about the rationale behind the program and its features. 

## How to cite

`MerquryRising` has not yet been published. Please cite this GitHub in the meantime.

## Installation requirements

1. Clone the GitHub repo.
2. Install RStudio with the appropriate packages, below. (RStudio will probably ask you if you want to install these when you open the `MerquryRising.Rmd` file.)
3. Open the `MerquryRising.Proj` in RStudio and then open the `MerquryRising.Rmd` file.
4. Test-knit the Rmarkdown.

`MerquryRising` uses the following libraries:

```
library(tidyverse)
library(ggridges)
library(grDevices)
library(GGally)
library(RColorBrewer)
library(ggstatsplot)
library(writexl)
library(kableExtra)
library(tools)
```

## How to run

1. Setup a directory with the relevant `merqury` outputs. (See GitHub repo for an example.)
2. Copy the `MerquryRising.Rmd` and update the settings as appropriate (or edit the config file).
3. If needed/wanted, create a `merquryrising.fofn` file to establish aliases for the assemblies and set the ordering to be used by `histsort`.
4. Knit the Rmd.
5. Optionally, edit the Rmd and re-knit until you are happy.

The Merqury files needed are:
* Read kmer `*.hist.ploidy` file.
* A set of `*.only.hist` files of assembly-only kmer counts.
* A set of `*.spectra-cn.hist` files of read and assembly kmer counts.

## Troubleshooting

If you are not getting the outputs you expect, please check or update the `merquryrising.config` configuration file, consisting of tab-delimited pairs of `setting` and `value`. If things are not working, please raise an issue on GitHub!