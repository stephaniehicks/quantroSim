# quantroSim

The R-package **quantroSim** is the supporting 
   data simulation package for the **quantro** 
   R/Bioconductor package which can be used to 
   simulate data from the epigenome including
   gene expression (microarrays, but will be extended 
   to include RNA-Seq) and DNA methylation (microarrays).  
For help with the **quantroSim** R-package, there is a vignette
  available in the /vignettes folder.
  
### Installation

 To install the package, you can check out this Github repository and install from source or use the `install_github()` function in the **devtools** R package:
```s
library(devtools)
install_github(repo = "quantroSim", username = "stephaniehicks")
```

After installation, the package can be loaded into R.
```s
library(quantroSim)
```

