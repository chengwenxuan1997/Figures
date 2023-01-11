# Preparation -------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
work.path <- ifelse(test = is.na(args[1]),
                    yes = "d:/genomicdata/External/Xlu/TimingLandmark",
                    no = args[1])
Tumor.file <- ifelse(test = is.na(args[2]),
                     yes = "d:/genomicdata/External/Xlu/CopyNumber/InputData/Bam/KAP2T.paired.chr6.bam",
                     no = args[2])
Normal.file <- ifelse(test = is.na(args[3]),
                      yes = "d:/genomicdata/External/Xlu/CopyNumber/InputData/Bam/PBL2.paired.chr6.bam",
                      no = args[3])
genome <- ifelse(test = is.na(args[4]),
                 yes = "hg38",
                 no = args[4])
# work.path <- "I:/genomicdata/External/Xlu/TimingLandmark/";setwd(work.path)
article.path <- file.path(work.path, "Literatures")
data.path <- file.path(work.path, "InputData")
pkg.path <- file.path(work.path, "Packages")
code.path <- file.path(work.path, "Codes")
res.path <- file.path(work.path, "Results")
fig.path <- file.path(work.path, "Figures")
ref.path <- file.path(work.path, "Reference")

invisible(lapply(ls()[grep("path", ls())], function(x) if (!dir.exists(get(x))) dir.create(get(x))))

cat(work.path, "\n")

# BiocManager::install("igordot/copynumber")
# devtools::install_github("Wedge-Oxford/battenberg")

library(Battenberg)

Reference <- read.table(file.path(data.path, paste0(genome, ".ReferenceList.txt")), header = T)
Reference <- setNames(object = file.path(ref.path, genome, Reference$Filename), nm = Reference$Name)

filename <- file.path(ref.path, genome, "impute_info.txt")
file.copy(from = paste0(filename, ".bak"), to = filename, overwrite = T)
xfun::gsub_file(file = filename, pattern = "REF_PATH", replacement = file.path(ref.path, genome))

loci.prefix <- switch(genome,
                      "hg38" = "1kg.phase3.v5a_GRCh38nounref_loci_chr",
                      "hg19" = "1000genomesloci2012_chr")
allele.prefix <- switch(genome,
                        "hg38" = "1kg.phase3.v5a_GRCh38nounref_allele_index_chr",
                        "hg19" = "1000genomesAlleles2012_chr")

args <- list(analysis = "paired",
             tumourname = "Tumor", normalname = "Normal",
             tumour_data_file = Tumor.file,
             normal_data_file = Normal.file,
             imputeinfofile = file.path(ref.path, genome, "impute_info.txt"),
             g1000prefix = file.path(Reference["1000G_loci"], loci.prefix),
             g1000allelesprefix = file.path(Reference["1000G_loci"], allele.prefix),
             problemloci = file.path(Reference["probloci"]),
             gccorrectprefix = file.path(Reference["GC_correction"], "1000G_GC_chr"),
             repliccorrectprefix = file.path(Reference["RT_correction"], "1000G_GC_chr"),
             ismale = T,
             data_type = "wgs",
             impute_exe = "impute2",
             allelecounter_exe = "alleleCounter",
             nthreads = 4,
             platform_gamma = 1,
             phasing_gamma = 1,
             segmentation_gamma = 10,
             segmentation_kmin = 3,
             phasing_kmin = 1,
             clonality_dist_metric = 0,
             ascat_dist_metric = 1,
             min_ploidy = 1.6, max_ploidy = 4.8,
             min_rho = 0.4, min_goodness = 0.63,
             uninformative_BAF_threshold = 0.51,
             min_normal_depth = 10, min_base_qual = 20, min_map_qual = 35,
             calc_seg_baf_option = 3,
             skip_allele_counting = F, skip_preprocessing = F, skip_phasing = F,
             externalhaplotypefile = NA,
             usebeagle = F, beaglejar = NA, beagleref.template = NA, beagleplink.template = NA,
             beaglemaxmem = 10, beaglenthreads = 1, beaglewindow = 40, beagleoverlap = 4,
             javajre = "java", write_battenberg_phasing = T,
             multisample_relative_weight_balanced = 0.25, multisample_maxlag = 100,
             segmentation_gamma_multisample = 5,
             snp6_reference_info_file = NA,
             heterozygousFilter = "none", prior_breakpoints_file = NULL,
             GENOMEBUILD = genome)
args <- args[intersect(names(args), formalArgs(battenberg))]
do.call(battenberg, args)
