# The code for the manifest creation is based on the code included in the package
# IlluminaHumanMethylationEPICanno.ilm10b4.hg19.
# Raw files accessed 12 Aug 2022 from the Illumina Website (Package documentation)
# MD5 (204617710004_R01C01_Grn.idat) = ccd965755e03c13dd846b46f516d6e7b
# MD5 (MouseMethylation-12v1-0_A2_fixed_probetype.csv) = c1e1c06d1564b031a86cf14e8623296d

library(minfi)
library(illuminaio)
library(devtools)
library(dplyr)
library(here)
library(tidyverse)
source("inst/scripts/read.manifest.R") # adapted read manifest file from minfi; there are a few duplicated probes (7786/2). Use Illumina IDs rather than names here.

idat.filepath <- "inst/data/204617710004_R01C01_Grn.idat"
manifest.filepath <- "inst/data/MouseMethylation-12v1-0_A2_fixed_probetype.csv"

if(!file.exists(idat.filepath) || !file.exists(manifest.filepath)) {
  cat("Missing files, quitting\n")
  q(save = "no")
}

maniTmp <- read.manifest(manifest.filepath)
anno <- maniTmp$manifest
manifestList <- maniTmp$manifestList

## Checking
mouse <- readIDAT(idat.filepath)
address.mouse <- as.character(mouse$MidBlock)
dropCpGs <- anno$Name[anno$AddressB_ID != "" & !anno$AddressB_ID %in% address.mouse]
dropCpGs <- anno$Name[anno$AddressA_ID != "" & !anno$AddressA_ID %in% address.mouse]
table(substr(dropCpGs, 1,2))

## probeInfo for ChAMP package
probeInfo <- anno %>%
  select(Name, Infinium_Design_Type) %>%
  mutate(Infinium_Design_Type = case_when(Infinium_Design_Type == "I" ~ 1,
                                          Infinium_Design_Type == "II" ~ 2))
## Manifest package
IlluminaMouseMethylationmanifest <- do.call(IlluminaMethylationManifest,
                                            list(TypeI = manifestList$TypeI,
                                                 TypeII = manifestList$TypeII,
                                                 TypeControl = manifestList$TypeControl,
                                                 TypeSnpI = manifestList$TypeSnpI,
                                                 TypeSnpII = manifestList$TypeSnpII,
                                                 annotation = "IlluminaMouseMethylation"))


stopifnot(validObject(IlluminaMouseMethylationmanifest))

save(IlluminaMouseMethylationmanifest, compress = "xz",
     file = "data/IlluminaMouseMethylationmanifest.rda")

save(probeInfo, compress = "xz",
     file = "data/probeInfoMouse.rda")
