# scShaper

A tool for linear trajectory inference from **s**ingle-**c**ell RNA-seq data through discrete pseudotime estimation using **S**hortest **HA**miltonian path **PER**muted clustering.

## Description

scShaper is an R package implementing a new trajectory inference method, which enables accurate and fast linear trajectory inference from single-cell RNA-seq (scRNA-seq) data. The method is based on estimating discrete pseudotimes using shortest Hamiltonian path permuted clusterings. The clustering is performed using the *k*-means algorithm, and the result is permuted using a special case of Kruskal's algorithm. scShaper uses Principal Component Analysis (PCA) to analyze linear dependencies between the pseudotimes. Based on the PCA results, scShaper selects the pseudotimes that contribute most to the first princal component and builds an ensemble pseudotime based on these by flipping and averaging the selected pseudotimes. Finally, scShaper smooths the ensemble pseudotime using local regression (LOESS).

## Installation

The easiest way to install scShaper is to use devtools R package.

```R
devtools::install_github("elolab/scshaper")
```

## Examples

[The first example](https://htmlpreview.github.io/?https://github.com/elolab/scshaper/blob/main/examples/1_scShaper_with_dyno_integration.nb.html) demonstrates scShaper using a simulated scRNA-seq dataset and integration with the dyno pipeline. This workflow is recommended if the user wants to use scShaper with **scRNA-seq data**. 

[The second example](https://htmlpreview.github.io/?https://github.com/elolab/scshaper/blob/main/examples/2_scShaper_mathematical_trajectory.nb.html) demonstrates scShaper with a simulated trigonometric dataset (2D spiral). This workflow is recommended if the user wants to use scShaper with **other data types that require no dimensionality reduction**. 

## Citation

*Johannes Smolander. Sini Junttila. Mikko S. Venäläinen. Laura L. Elo. scShaper: ensemble method for fast and accurate linear trajectory inference from single-cell RNA-seq data. GitHub. [https://github.com/elolab/scshaper.](https://github.com/elolab/scshaper.)*
