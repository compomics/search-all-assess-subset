# search-all-assess-subset
An implementation of the Search All, Asses Subset strategy for FDR estimation in shotgun proteomics.

<!-- A preprint of the manuscript on the search-all-assess FDR procedure can be found on bioRxiv ([https://doi.org/10.1101/094581](https://doi.org/10.1101/094581)). -->
The search-all-assess-subset FDR procedure has been published in [Sticker et al. Nature Methods 14, 643–644 (2017) doi:10.1038/nmeth.4338](https://doi.org/10.1038/nmeth.4338)
>Mass spectrometrists should search for all peptides, but assess only the ones they care about<br/>
>(Adriaan Sticker, Lennart Martens, Lieven Clement)

Please include this reference if you use search-all-assess-subset FDR procedure as part of a publication.

A user-friendly GUI version of saas is hosted [here](http://iomics.ugent.be/saas/).

More information on the FDR estimation procedure for subsets and the usage of the SAAS package can be found in the [vignette](https://compomics.github.io/search-all-assess-subset/articles/saas.html)
<!-- (http://htmlpreview.github.io/?https://github.com/compomics/search-all-assess-subset/blob/master/vignettes/saas.html). -->

## Installation

 You first need to install the [devtools](https://cran.r-project.org/package=devtools) package.

```r
install.packages("devtools")
```

Then install saas using the `install_github` function in the
[devtools](https://cran.r-project.org/package=devtools) package.
```r
library(devtools)
install_github("compomics/search-all-assess-subset")
```

## Running the saas GUI

When saas is installed, you can easily launch the GUI from R.
By default, you can see the GUI by visiting [http://127.0.0.1:3320](http://127.0.0.1:3320)
```r
library(saas)
saas_gui()
```
Or if you want to run saas immediately from the terminal.
```r
 R -e "library(saas);saas_gui()"
```
