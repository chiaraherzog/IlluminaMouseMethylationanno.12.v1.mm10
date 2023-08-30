# The code for the manifest creation is based on the code included in the package
# IlluminaHumanMethylationEPICanno.ilm10b4.hg19.
# Raw files accessed 12 Aug 2022 from the Illumina Website (Package documentation), not uploaded to repository due to file size limitations

library(minfi)
library(illuminaio)
library(devtools)
library(dplyr)
library(here)
library(tidyverse)
source("inst/scripts/read.manifest.R") # adapted read manifest file from minfi; there are a few duplicated probes (7786/2). Use Illumina IDs rather than names here.

idat.filepath <- "inst/data/204617710004_R01C01_Grn.idat"
manifest.filepath <- "inst/data/MouseMethylation-12v1-0_A2.csv"

if(!file.exists(idat.filepath) || !file.exists(manifest.filepath)) {
  cat("Missing files, quitting\n")
  q(save = "no")
}

# Fix Illumina Probe types - no 2/1 allowed, only "II" or "I"
maniTmp <- read.table(manifest.filepath, sep = ",", skip = 7, header = T) %>%
  mutate(Infinium_Design_Type = case_when(Infinium_Design_Type == 2 ~ "II",
                                          Infinium_Design_Type == 1 ~ "I")) %>%
  write.table(file = "inst/data/MouseMethylation-12v1-0_A2_fixed_probetype.csv",
              sep = ",", quote = F, row.names = F)

manifest.filepath <- "inst/data/MouseMethylation-12v1-0_A2_fixed_probetype.csv"
maniTmp <- read.manifest(manifest.filepath)
anno <- maniTmp$manifest

# transform Islands 
Islands <- anno %>%
  tidyr::pivot_longer(c(N_Shelf, N_Shore, CpG_Island, S_Shore, S_Shelf),
                      names_to = "Relation_to_Island",
                      values_to = "Isl2") %>%
  group_by(Name, Relation_to_Island) %>%
  summarise(Island = paste0(Isl2)) %>%
  filter(Island != "") %>%
  group_by(Name) %>%
  mutate(Relation_to_Island = unique(Relation_to_Island)[1], # keep first anno
         Island = paste0(unique(Island), collapse = ";")) %>%
  dplyr::mutate(Relation_to_Island = gsub("CpG_Island", "Island", Relation_to_Island)) |> 
  distinct()

# add open seas
library(tidyverse)
anno <- Islands %>%
  dplyr::full_join(anno, by = "Name") %>%
  mutate(Relation_to_Island = case_when(is.na(Relation_to_Island) ~ "OpenSea",
                                        TRUE ~ Relation_to_Island))

manifestList <- maniTmp$manifestList

## Checking
mouse <- readIDAT(idat.filepath)
address.mouse <- as.character(mouse$MidBlock)
dropCpGs <- anno$Name[anno$AddressB_ID != "" & !anno$AddressB_ID %in% address.mouse]
dropCpGs <- anno$Name[anno$AddressA_ID != "" & !anno$AddressA_ID %in% address.mouse]
table(substr(dropCpGs, 1,2))


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


## Annotation package
geneAnno.file <- "inst/data/MouseMethylation-12v1-0_A1_Annotation_Mus_musculus.csv"
geneAnno <- read.table(geneAnno.file, sep = ",", header = T)

geneAnno.UCSC <- geneAnno %>%
  rename(cg = name) %>%
  filter(Source == "UCSC") %>%
  group_by(cg) %>%
  summarise(Gene.UCSC = paste(unique(Gene[!is.na(Gene)]), collapse = ";"),
            Transcript.UCSC = paste(unique(Transcript[!is.na(Transcript)]), collapse = ";"),
            Feature.UCSC = paste(unique(Feature[!is.na(Feature)]), collapse = ";"))

geneAnno.NCBI <- geneAnno %>%
  rename(cg = name) %>%
  filter(Source == "NCBI") %>%
  group_by(cg) %>%
  dplyr::summarise(Gene.NCBI = paste(unique(Gene[!is.na(Gene)]), collapse = ";"),
                   Transcript.NCBI = paste(unique(Transcript[!is.na(Transcript)]), collapse = ";"),
                   Feature.NCBI = paste(unique(Feature[!is.na(Feature)]), collapse = ";"))

anno <- anno %>%
  full_join(geneAnno.UCSC) %>%
  full_join(geneAnno.NCBI)

nam <- names(anno)
names(nam) <- nam
nam[c("AddressA_ID", "AddressB_ID", "AlleleA_ProbeSeq", "AlleleB_ProbeSeq",
      "Infinium_Design_Type", "Next_Base", "Color_Channel", "Rep_Num")] <-  c("AddressA", "AddressB",
                                                                              "ProbeSeqA", "ProbeSeqB",
                                                                              "Type", "NextBase", "Color", "RepNum")


names(nam) <- NULL
names(anno) <- nam
anno <- data.frame(anno)
rownames(anno)<- anno$Name
Locations <- anno[, c("CHR", "MAPINFO")]
names(Locations) <- c("chr", "pos")
Locations$pos <- as.integer(Locations$pos)
Locations$chr <- paste("chr", Locations$chr, sep = "")
Locations$strand <- ifelse(anno$Strand == "F", "+", "-")
table(Locations$chr, exclude = NULL)
Locations <- as(Locations, "DataFrame")

Manifest <- anno[, c("Name", "AddressA", "AddressB",
                     "ProbeSeqA", "ProbeSeqB", "Type", "NextBase", "Color", "RepNum")]
Manifest <- as(Manifest, "DataFrame")

Islands.UCSC <- anno[, c("Island", "Relation_to_Island")]
names(Islands.UCSC) <- c("Islands_Name", "Relation_to_Island")
Islands.UCSC <- as(Islands.UCSC, "DataFrame")
table(Islands.UCSC$Relation_to_Island, exclude = NULL)

GenesUCSC <- anno[, c("Gene.UCSC", "Feature.UCSC", "Transcript.UCSC")]
names(GenesUCSC) <- c("GeneName_UCSC", "Feature_UCSC", "Transcript_UCSC")
GenesUCSC <- as(GenesUCSC, "DataFrame")

GenesNCBI <- anno[, c("Gene.NCBI", "Feature.NCBI", "Transcript.NCBI")]
names(GenesNCBI) <- c("GeneName_NCBI", "Feature_NCBI", "Transcript_NCBI")
GenesNCBI <- as(GenesNCBI, "DataFrame")

annoNames <- c("Locations", "Manifest", "Islands.UCSC", "GenesUCSC", "GenesNCBI")

for(nam in annoNames) {
  cat(nam, "\n")
  save(list = nam, file = file.path("data/", paste(nam, "rda", sep = ".")), compress = "xz")
}

annoStr <- c(array = "IlluminaMouseMethylation",
             annotation = "12.v1",
             genomeBuild = "mm10")

defaults <- c("Locations", "Manifest", "Islands.UCSC", "GenesNCBI", "GenesUCSC")
pkgName <- sprintf("%sanno.%s.%s", annoStr["array"], annoStr["annotation"],
                   annoStr["genomeBuild"])

annoObj <- IlluminaMethylationAnnotation(objectNames = annoNames, annotation = annoStr,
                                         defaults = defaults, packageName = pkgName)

assign(pkgName, annoObj)
save(list = pkgName,
     file = file.path("data/", paste(pkgName, "rda", sep = ".")), compress = "xz")
sessionInfo()

