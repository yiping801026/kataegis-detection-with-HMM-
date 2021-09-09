This file includes the optimized HMM function 'performHMM_op()'
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
performHMM_op <- function(x, chr_num){
   
  ### import library
  library(ggplot2)
  library(depmixS4)
  library(KatDetect)
 
  ### settings for HMM
  minMutsInKat = 6
  nstates = 3
  ### import data + log change.
 
  # Filter on selected chromosome.
  subset.x <- x[seqnames(x) %in% chr_num]
 
  # Filter IMD
  x.IMD <- intermutDistance(x)
  x.IMD <- x.IMD[as.vector(seqnames(x) %in% chr_num)]
  subset.x@intermutDistance <- x.IMD
 
  # Filter mutContexts
  x.mutContext <- mutContext(x)
  x.mutContext <- x.mutContext[as.vector(seqnames(x) %in% chr_num)]
  subset.x@mutContext <- x.mutContext
 
  if (length(subset.x@intermutDistance) < 10 ){return(NULL)}
 
  # Convert to tibble.
  logtran_data <- tibble::tibble(
    sampleNames = as.character(sampleNames(subset.x)),
    chr = as.character(seqnames(subset.x)),
    start = GenomicRanges::start(subset.x),
    end = GenomicRanges::end(subset.x),
    IMD = subset.x@intermutDistance,
    logdata = log10(IMD),
    mutContext = subset.x@mutContext
  )
 
  raw_logtran_data <- logtran_data %>% tibble::rowid_to_column()
 
  # Remove NA values (last variant has no subsequent variant)
  logtran_data <- logtran_data %>% dplyr::filter(!is.na(logdata))
 
  # Create lookup-index
  logtran_data <- logtran_data %>% tibble::rowid_to_column()
 
  # Get APOBEC mutational context.
  logtran_data <- logtran_data %>% dplyr::mutate(
    mutContext.APO = substring(mutContext,3,5),
    isAPO = ifelse(mutContext.APO %in% c('C>T', 'C>G'), T, F)
  )
 
  ### fit model using HMM
  ### change the data into categorical data '1' is low, '2' is middle and '3' is high.
  logtran_data <- logtran_data %>%
    dplyr::mutate(group = dplyr::case_when(
      IMD <= 3000 ~ 1,
      IMD > 3000 & IMD < 10000 ~ 2,
      IMD >= 10000 ~ 3))
 
  # Check if at least two groups are present, else return NULL.
  if(dplyr::n_distinct(logtran_data$group) < 2) return(NULL)
 
  # Check if all three states are present.
  if(dplyr::n_distinct(logtran_data$group) == 3 ){
   
    mod <- depmixS4::depmix(
      group ~ 1,
      data = logtran_data,
      nstates = nstates,
      family = multinomial('identity'),
      instart  = c(1, 0, 0),
      trstart = c(0.995, 0.001, 0.004, 0.01, 0.98, 0.01, 0.01, 0.01, 0.98),
      respstart = c(0.04, 0.03, 0.93, 0.12, 0.20, 0.68, 0.99, 0.01, 0.0)
    )
  }
 
  # If two states and two groups are present.
  if(dplyr::n_distinct(logtran_data$group) != 3){
   
    mod <- depmixS4::depmix(
      group ~ 1,
      data = logtran_data,
      nstates = 2,
      family = multinomial('identity'),
      instart  = c(1, 0),
      trstart = c(0.995, 0.005, 0.02, 0.98),
      respstart = c(0.01, 0.99, 1.00, 0.00)
    )
  }
 
  fit.mod <- depmixS4::fit(mod,emcontrol = depmixS4::em.control(maxit = 1, random.start = F), verbose = F)
  summary(fit.mod)
 
  ### set.states is a table includes the results of HMM  
  set.states <- depmixS4::posterior(fit.mod)
  logtran_data$state <- set.states$state
 
  ### make a table to distinguish whether 1.the chromosome has kataegis state 2.the state is kataegis or not    
  ### this table includes details of every continuous region/segment of one state along the chromosome
  ### 1. 'start' mutation index; 2. 'end' mutation index; 3. 'state' of this region; 4. 'mean_IMD' means the average distance of the dots in this region; 5. 'K' means the number of dots in this region; 6. 'p_APO' means the percentage of mutations from APOBEC enzymes among this region. 7. 'kat' means the states is kataegis state or not.
 
  logtran_data <- logtran_data %>% dplyr::mutate(segment = state != c(state[-1], NA))
 
  # Give each segment a unique identifier.
  logtran_data <- logtran_data %>% mutate(ticker = cumsum(segment))
  logtran_data$segmentId <- c(0, logtran_data$ticker[1:nrow(logtran_data)-1])
  logtran_data$ticker <- NULL
 
  # State-wise statistics.
  logtran_data <- logtran_data %>%
    dplyr::group_by(state) %>%
    dplyr::mutate(mean_IMD_State = mean(IMD))
 
  # Segment-wise statistics.
  logtran_data <- logtran_data %>%
    dplyr::group_by(segmentId) %>%
    dplyr::mutate(
      nMuts = dplyr::n()+1,
      mean_IMD = mean(IMD),
      num_APO = sum(isAPO),
      p_APO = num_APO / dplyr::n(),
      startMutOfSegment = min(rowid),
      endMutOfSegment = max(rowid),
      mean_IMD_State = unique(mean_IMD_State),
      state = unique(state)
    ) %>%
    dplyr::ungroup()
 
  ### Whether the chromosome have kataegis state or not. If 'mm_IMD ', the mean of mean_IMD for all regions of one state smaller than 3000bp
  ### For the chromosome includes kataegis state, distinguish which state is kataegis state
  logtran_data <- logtran_data %>%
    dplyr::group_by(state) %>%
    dplyr::mutate(
      isKataegisState = ifelse(any(mean_IMD_State <= 3000 & nMuts >= minMutsInKat), T, F),
      isKataegisSegment = ifelse(mean_IMD <= 3000 & nMuts >= minMutsInKat, T, F),
    ) %>%
    dplyr::ungroup()
 
  if(any(logtran_data$isKataegisState)) print(sprintf('There is kataegis state in this chromosome! State: %s', unique(logtran_data %>% dplyr::filter(isKataegisState) %>% dplyr::pull(state))))
 
  if(any(logtran_data$isKataegisSegment)) print(sprintf('HMM found %s kataegis segment(s) in this chromosome!', dplyr::n_distinct(logtran_data %>% dplyr::filter(isKataegisSegment) %>% dplyr::pull(segmentId))))
 
  # Print the kataegis segments.
  logtran_data %>% dplyr::filter(isKataegisSegment) %>% dplyr::distinct(segmentId, startMutOfSegment, endMutOfSegment, mean_IMD, nMuts, num_APO, p_APO)
 
  # change the order of chromosomes by 1:22 + X
  logtran_data$chr = factor(logtran_data$chr, paste0('chr',c(1:22,'X')),ordered = T)
  logtran_data = logtran_data[with(logtran_data,order(chr)),]
 
 
    p <- logtran_data %>%
    dplyr::mutate(stateClean = ifelse(!isKataegisState, sprintf('%s (Other)', state), sprintf('%s (Kataegis State)', state))) %>%
    ggplot2::ggplot(., aes(x = rowid, y = IMD, color = stateClean, alpha = isKataegisSegment, shape = isKataegisSegment)) +
    ggplot2::geom_point() +
    ggplot2::scale_shape_manual(values = c(1, 16)) +
    ggplot2::scale_alpha_manual(values = c(.5, 1)) +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000,100000000), limits = c(0, 100000000), labels = c('0bp', '10bp', '100bp', '1kb', '10kb', '100kb', '1Mb', '10Mb','100Mb'), expand = c(0,0)) +
    ggplot2::labs(x = 'Somatic mutations', y = 'Intermutational distance ') +
    ggplot2::theme(
      text = ggplot2::element_text(size=9, family='Helvetica', face = 'bold'),
      axis.text.x = ggtext::element_markdown(),
      axis.title.x = ggtext::element_textbox_simple(width = NULL, halign = .5),
      axis.title.y = ggtext::element_textbox_simple(size = 8, orientation = 'left-rotated', width = NULL, halign = .5),
      strip.text = ggtext::element_textbox_simple(width = NULL, halign = .5),
      panel.grid.major.x = ggplot2::element_line(colour = 'grey90', linetype = 'dotted'),
      panel.grid.major.y = ggplot2::element_line(colour = '#E5E5E5', linetype = 'dotted'),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
      panel.border = ggplot2::element_rect(fill = NA, colour = NA),
      strip.background = ggplot2::element_rect(colour = 'black', fill = 'white'),
      legend.text = ggtext::element_markdown()
    ) +
    ggplot2::facet_grid(.~chr, scales = 'free_x', space = 'free_x')
 
  print(p)
 
  ### return a table 'katMuts_Table' including sampleNames, chromosome, start and end location, mean IMD value for the segment, percentage of APOBEC, number of mutations in the segment.
 
  # TM is the transition matrix of HMM_op
  TM <- t(matrix(fit.mod@trDens,nrow = nstates,ncol = nstates))
  rownames(TM) <- c('state1', 'state2', 'state3')
  colnames(TM) <- c('state1', 'state2', 'state3')
 
  raw_logtran_data <- full_join(raw_logtran_data,logtran_data)
 
  logtran_data_IMD <- logtran_data %>%
    filter(isKataegisState == TRUE) %>%
    filter(isKataegisSegment  == TRUE)
 
  kat_segmentId <- unique(logtran_data_IMD$segmentId)
 
  list_katTable <- lapply(kat_segmentId, function(x){
    l <- logtran_data_IMD[logtran_data_IMD$segmentId == x,]$rowid
    katMuts <- tibble::tibble(
      sampleNames = as.character(logtran_data[logtran_data$rowid == l[1],]$sampleNames),
      chr = as.character(logtran_data[logtran_data$rowid == l[1],]$chr),
      start = logtran_data[logtran_data$rowid == l[1],]$end,
      end = raw_logtran_data[raw_logtran_data$rowid == l[length(l)]+1,]$start,
      mean_IMD = logtran_data[logtran_data$rowid == l[1],]$mean_IMD,
      p_APO = logtran_data[logtran_data$rowid == l[1],]$p_APO,
      nMuts = logtran_data[logtran_data$rowid == l[1],]$nMuts
    )
    return(katMuts)
  })
 
  katMuts_Table <- do.call(rbind, list_katTable)
 
  return(katMuts_Table)
 
}
```


From: Y. Ping
Sent: Tuesday, June 29, 2021 7:53 PM
To: Yi Ping <pingyi900@gmail.com>
Subject: codes_HMM_CPCT.R
 
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
