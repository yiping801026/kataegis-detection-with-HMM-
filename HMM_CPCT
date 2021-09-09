#use performHMM_op() function to check all samples of breast cancer on CPCT-02-DRUP dataset

### import library
library(KatDetect)

### import dataset--------------------------------------------------------
metadata.DR41 <- readr::read_tsv('/mnt/data2/hartwig/DR41/Feb2021/dataHMF/metadata/metadata.tsv')
Tbreast_samples <- metadata.DR41[which(metadata.DR41$primaryTumorLocation == 'Breast'),]

genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

# Retrieve VCF files.
path_table <- tibble::tibble(
    sampleId = Tbreast_samples$sampleId,
    sample = sprintf('KatRanges.%s', Tbreast_samples$sampleId),
    pathVCF = sprintf('/mnt/data2/hartwig/DR41/Feb2021/dataHMF/combinedData/%s.purple.somatic.vcf.gz', Tbreast_samples$sampleId)
)

### Import VCF files
samples_list <- BiocParallel::bplapply(path_table$pathVCF, function(x){
    KatDetect::importVCFasKatRanges(pathVCF = x, genome, field.AF = 'PURPLE_AF', keepSamples = gsub('\\..*', '', basename(x)), passOnly = T, keepAnnotation = NA , snpOnly = T)
}, BPPARAM = BiocParallel::MulticoreParam(workers = 1))

# Add names.
names(samples_list) <- unlist(lapply(samples_list, function(x) unique(Biobase::sampleNames(x))))

### perform optimized HMM on the dataset-------------------------------

set.seed(2000)

### Perform HMM per chromosome

results <- BiocParallel::bplapply(samples_list, function(x){

    # Perform HMM per chromosome.
    hmm.PerChr <- lapply(GenomeInfoDb::seqlevelsInUse(x), function(y){
        ParallelLogger::logInfo(sprintf('Working on: %s - %s', levels(sampleNames(x)), y))
        return(performHMM_op(x, chr_num = y))
    })
    #names(hmm.PerChr) <- GenomeInfoDb::seqlevelsInUse(x)
}, BPPARAM = BiocParallel::MulticoreParam(workers = 1, progressbar = T))
