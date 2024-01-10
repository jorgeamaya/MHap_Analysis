
if(!require(adegenet)){
  install.packages("adegenet")
  library(adegenet)
}

if(!require(ade4)){
  install.packages("ade4")
  library(ade4)
}

if(!require(poppr)){
  install.packages("poppr")
  library(poppr)
}

if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}

if(!require(magrittr)){
  install.packages("magrittr")
  library(magrittr)
}

if(!require(tidyr)){
  install.packages("tidyr")
  library(tidyr)
}

if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}

if(!require(cowplot)){
  install.packages("cowplot")
  library(cowplot)
}


if(!require(vegan)){
  install.packages("vegan")
  library(vegan)
}

if(!require(parallel)){
  install.packages("parallel")
  library(parallel)
}

if(!require(ape)){
  install.packages("ape")
  library(ape)
}

if(!require(pegas)){
  install.packages("pegasn")
  library(pegas)
}

if(!require(RColorBrewer)){
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

if(!require(Hmisc)){
  install.packages('Hmisc')
  library(Hmisc)
}

if(!require(ggpubr)){
  install.packages('ggpubr')
  library(ggpubr)
}

if(!require(doMC)){
  install.packages('doMC')
  library(doMC)
}

if(!require(svMisc)){
  install.packages('svMisc')
  library(svMisc)
}

if(!require(Biostrings)){
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("Biostrings")
  library(Biostrings)
}

if (!require("argparse", quietly = TRUE)) {
  install.packages("argparse")
  library(argparse)
}
