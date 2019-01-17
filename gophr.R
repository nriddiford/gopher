# if (!requireNamespace("BiocManager"))
#   install.packages("BiocManager")
# BiocManager::install("goseq", version = "3.8")
# BiocManager::install("org.Dm.eg.db", version = "3.8")

library(goseq)
library(plyr)
library(dplyr)
library(stringr)
library(data.table)
library(scales)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(pheatmap)
library(GO.db)
library(forcats)
library(GSEABase)

cleanTheme <- function(base_size = 12) {
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.spacing = unit(2, "lines"),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 30),
    plot.margin = unit(1:4, "line")
  )
}


get_data <- function(..., in_file = NULL, mut_type='sv'){
  if(mut_type=='sv'){
    cat("Reading SV data")
    if(missing(in_file))
      in_file = '~/Desktop/parserTest/filtered_231018/summary/merged/all_genes_filtered.txt'
    genes <- read.delim(in_file)
    colnames(genes) <- c("event", "sample", "genotype", "type", "af", "chrom", "gene")
    # genes <- genes %>% 
    #   dplyr::filter(...)
    
  } else if(mut_type=='snv'){
    cat("Reading SNV data")
    if(missing(in_file))
      in_file = '~/Desktop/script_test/mutationProfiles/data/annotated_snvs.txt'
    genes <- read.delim(in_file)
    genes <- genes %>% 
      # dplyr::filter(...) %>% 
      dplyr::mutate(chrom = chromosome, 
                    af = allele_frequency) %>% 
      dplyr::select(sample, chrom, af, status, gene)
    # colnames(genes) <- c("event", "sample", "genotype", "type", "af", "chrom", "gene")
  } else if(mut_type=='indel'){
    cat("Reading indel data")
    if(missing(in_file))
      in_file = '~/Desktop/script_test/mutationProfiles/data/annotated_indels.txt'
    genes <- read.delim(in_file)
    genes <- genes %>% 
      # dplyr::filter(...) %>% 
      dplyr::mutate(chrom = chromosome, 
                    af = allele_frequency) %>% 
      dplyr::select(sample, chrom, af, gene, status, type)
    # colnames(genes) <- c("event", "sample", "genotype", "type", "af", "chrom", "gene")
  }
  
  genes <- genes %>% 
    dplyr::mutate(cell_fraction = ifelse(chrom %in% c('X', 'Y'), af,
                                         ifelse(af*2>1, 1, af*2))) %>%
    dplyr::filter(...) %>% 
    dplyr::distinct(sample, gene, .keep_all=TRUE) %>% 
    droplevels()
  
  return(genes)
}

go_getter <- function(..., gene_uni='~/Desktop/script_test/svBreaks/inst/extdata/gene_lengths.txt',
                     hit_genes=NULL,
                     excluded_genes=NULL,
                     go_class=c("GO:BP", "GO:MF", "GO:CC"), by_sample=TRUE, p_threshold=0.005, mut_type='sv') {
  
  if(missing(hit_genes)) stop("Option 'hit_genes' required")
  if(missing(excluded_genes)){
    excluded_genes <- c()
  }
  
  gene_lengths <- read.delim(gene_uni, header = T)
  
  geneUniverse <- gene_lengths %>%
    dplyr::filter(!gene %in% excluded_genes) %>% 
    dplyr::select(gene) %>% 
    droplevels()
  
  genes_per_sample <- function(filter_sample=NULL, genes=hit_genes, term_filt=go_class) {
    if(missing(filter_sample)){
      filter_sample <- levels(genes$sample)
      sample_name <- 'all'
    } else{
      sample_name <- filter_sample
    }
    
    genes <- genes %>%
      dplyr::filter(sample %in% filter_sample) %>%
      dplyr::distinct(gene) %>%
      dplyr::mutate(val = 1) %>%
      dplyr::select(gene, val)
    
    cat("Found", nrow(genes), "hit genes:", "\n", paste(genes$gene), "\n")
    nrow(geneUniverse)
    universe <- dplyr::anti_join(geneUniverse, genes)
    nrow(universe)
    
    universe$val <- 0
    all_genes <- rbind(genes, universe)
    
    gene_vals <- all_genes$val
    names(gene_vals) <- all_genes$gene
    
    pwf <- goseq::nullp(gene_vals,"dm3","geneSymbol", plot.fit = F)
    GO_wall <- goseq::goseq(pwf,"dm3","geneSymbol", test.cats = term_filt)
    
    over_rep_hits <- GO_wall$category[p.adjust(GO_wall$over_represented_pvalue, method="BH")<.05]
    under_rep_hits <- GO_wall$category[p.adjust(GO_wall$under_represented_pvalue, method="BH")<.05]
    
    cat(length(over_rep_hits), "significantly enriched terms", "\n")
    cat(length(under_rep_hits), "significantly depleted terms", "\n")
    
    filtered_go <- GO_wall %>% 
      dplyr::filter(over_represented_pvalue <= p_threshold || under_represented_pvalue <= p_threshold,
                    !is.na(term)) %>% 
      dplyr::select(over_represented_pvalue, under_represented_pvalue, term, category, ontology)
    
    cat(nrow(filtered_go), term_filt, "terms with p_val <=",  p_threshold, "\n")
    colnames(filtered_go) <- c(paste0(sample_name, "_over"), paste0(sample_name, "_under"), 'term', 'GO', 'ontology')
    
    return(list(filtered_go, over_rep_hits, under_rep_hits))
  }
  
  if(by_sample){
    cat("Running on each sample individually\n")
    samples <- levels(hit_genes$sample)
    
    sample_list <- list()
    
    for(s in samples){
      cat("GSEA on sample:", s, "\n")
      sample_list[[s]] <- genes_per_sample(filter_sample = s, genes=hit_genes)[[1]]
    }
    all_hits <- rbindlist(sample_list, fill = TRUE)
  } else{
    cat("Running on all samples together\n")
    all_hits <- genes_per_sample(genes=hit_genes)
  }
  return(all_hits)
}


go_figure <- function(..., combined=NULL, go_class=c("BP", "MF", "CC"), p_threshold = 1, n_samples = 10, direction = 'over', print=FALSE){
  if(missing(combined)) stop("Must provide dataframe with goSeq results. Exiting")
  if(!direction %in% c('over', 'under')) stop("Direction must be set to 'over' or 'under'. Exciting")
  cat("Plotting ", direction, "-represented terms\n", sep = '')
  replacer = 'under'
  if (direction == 'under'){
    replacer = 'over'
  }
  
  cleaned <- combined %>% 
    dplyr::filter(ontology %in% go_class) %>% 
    dplyr::select(-contains(replacer)) %>% 
    dplyr::mutate(term = gsub("'", '', term)) %>% 
    # dplyr::mutate(term = substr(term, start = 1, stop = 50)) %>% 
    dplyr::select(-ontology, -GO) %>% 
    droplevels()
  
  colnames(cleaned) <- gsub(x = names(cleaned), pattern = "_.*", replacement = "")
  
  df <- reshape2::melt(cleaned, id='term') %>% 
    dplyr::mutate(sample = as.factor(variable),
                  go_term = as.factor(term)) %>% 
    dplyr::select(-term, -variable)
  
  
  ### This rescales within go_term
  tallied_terms <- df %>%
    dplyr::filter(...) %>% 
    dplyr::mutate(value = ifelse(is.na(value), 1, value)) %>% 
    dplyr::group_by(sample) %>% 
    dplyr::arrange(value) %>% 
    dplyr::distinct(go_term, .keep_all=TRUE) %>% 
    
    # dplyr::filter(!is.na(value)) %>% 
    # dplyr::mutate(value = ifelse(is.na(value), 1, value)) %>% 
    dplyr::group_by(go_term) %>% 
  
    dplyr::mutate(count = sum(value < p_threshold)) %>% 
    # dplyr::filter(count >= n_samples) %>% 
    # dplyr::ungroup() %>% # To perform global rescale
    # dplyr::mutate(rescale = ifelse(count >= n_samples, rescale(-value), 0)) %>% 
    dplyr::filter(count >= n_samples) %>% 
    dplyr::ungroup() %>% # To perform global rescale
    dplyr::mutate(rescale = rescale(-value)) %>% 
    dplyr::arrange(-count)
  
  wide_tally <- tallied_terms %>%
    dplyr::select(go_term, sample, rescale) %>% 
    # unite(var, go_term, sample) %>% 
    tidyr::spread(sample, rescale)
  
  
  wt_copy <- wide_tally
  wide_tally$go_term <- NULL

  m1 <- as.matrix(wide_tally)
  
  rownames(m1) <- wt_copy$go_term

  if(direction=='under'){
    colours = rev(c("#A81B1B", "#CF5555", "#E39191", "#EDD0D0", "#F5EBEB"))
  } else{
    colours = c("#EBF4F5", "#7BB3D4", "#65ABC2", "#4594B3", "#366C9E")
  }
  pheatmap::pheatmap(m1, color = colours, cutree_rows = 3, cutree_cols = 3)
  
  if(print){
    return(tallied_terms)
  }
  
    # tidyr::spread(go_term, sample)
  
  # colnames(wide_tally) <- gsub(x = names(wide_tally), pattern = "rescale.", replacement = "")  
  
  
  # wide_tally <- reshape(tallied_terms, idvar = "go_term", timevar = "sample", direction = "wide")
  #   
  # 
  # # This is maybe better: 
  # #  Keep terms where at least one sample has a value <= p_threshold
  # ontologies_tallied <- df %>% 
  #   dplyr::group_by(go_term) %>% 
  #   dplyr::filter(!is.na(value),
  #                 any(value <= p_threshold) ) %>% 
  #   droplevels() %>% 
  #   dplyr::add_tally()
  # 
  # df2 <- plyr::ddply(ontologies_tallied, .(sample), transform,
  #                    rescale = rescale(-value))
  # 
  # df2 <- df2 %>% 
  #   # dplyr::filter(value <= p_threshold) %>% 
  #   dplyr::filter(n >= n_hits) %>% 
  #   dplyr::select(-n, -value) %>% 
  #   dplyr::arrange(-rescale, sample, go_term) %>% 
  #   droplevels()
  # 
  # df4 <- reshape(df2, idvar = "go_term", timevar = "sample", direction = "wide")
  # colnames(df4) <- gsub(x = names(df4), pattern = "rescale.", replacement = "")  
  # 
  # df5 <- df4
  # df5$term <- NULL
  # m1 <- as.matrix(df5)
  # rownames(m1) <- df4$term
  # m1[is.na(m1)] <- 0
  # 
  # pheatmap::pheatmap(m1, color = c("#EBF4F5", "#7BB3D4", "#65ABC2", "#4594B3", "#366C9E"), cutree_rows = 3, cutree_cols = 3)
  # 
  save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
    png(filename, width = width, height = height, res = res)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  # save_pheatmap_png(my_heatmap, "pheatmap_all_p_0.001.png")
}



gopher_broke <- function(..., combined, slim_file = 'goslim.obo', go_class='BP', direction='over', print=FALSE, evidenceCode="TAS"){
  if(!direction %in% c('over', 'under')) stop("Direction must be set to 'over' or 'under'. Exciting")
  
  replacer = 'under'
  if (direction == 'under'){
    replacer = 'over'
  }
  
  cleaned <- combined %>% 
    dplyr::filter(ontology == go_class) %>% 
    dplyr::select(-contains(replacer)) %>% 
    dplyr::mutate(term = gsub("'", '', term)) %>% 
    dplyr::select(-ontology, -term) %>% 
    droplevels()
  
  colnames(cleaned) <- gsub(x = names(cleaned), pattern = "_.*", replacement = "")
  
  df <- reshape2::melt(cleaned, id='GO') %>% 
    dplyr::mutate(sample = as.factor(variable),
                  go_term = GO) %>% 
    dplyr::select(-GO, -variable)
  
  filtered_terms <- df %>%
    # dplyr::filter(...) %>% 
    dplyr::mutate(value = ifelse(is.na(value), 1, value)) %>% 
    dplyr::group_by(sample) %>% 
    dplyr::arrange(value) %>% 
    dplyr::distinct(go_term, .keep_all=TRUE) %>% 
    dplyr::filter(value < 1)
  
  slims_per_sample <- list()
  for(s in levels(filtered_terms$sample)){
    cat(s, "\n")
    sample_terms <- filtered_terms %>%
      dplyr::filter(sample == s)
    
    go_terms <- sample_terms$go_term
    myCollection <- GOCollection(go_terms)
    
    slim <- getOBOCollection(slim_file)
    slim_counts <- goSlim(myCollection, slim, go_class, evidenceCode %in% evidenceCode)
    
    slims_per_sample[[s]] <- data.frame(sample = s, term=slim_counts$Term, percent = slim_counts$Percent, count=slim_counts$Count)
  }
  
  all_slims <- do.call(rbind, slims_per_sample)
  rownames(all_slims) <- NULL
  
  present_terms <- all_slims %>% 
    dplyr::filter(count > 0) %>% 
    dplyr::group_by(term) %>% 
    dplyr::mutate(across_samples = n()) %>% 
    dplyr::group_by(sample) %>%
    dplyr::mutate(score = rescale(percent)) %>% 
    dplyr::top_n(10, percent)
  
  wide_tally <- present_terms %>%
    dplyr::select(term, sample, percent) %>% 
    # unite(var, go_term, sample) %>% 
    tidyr::spread(sample, percent)
  
  wt_copy <- wide_tally
  wide_tally$term <- NULL
  wide_tally[is.na(wide_tally)] <- 0
  m1 <- as.matrix(wide_tally)
  
  rownames(m1) <- wt_copy$term
  
  breaksList = seq(0, 100, by = 20)
  
  if(direction=='under'){
    colours <- colorRampPalette(brewer.pal(n = 5, name = "Reds"))(length(breaksList))
  } else{
    colours <- colorRampPalette(brewer.pal(n = 5, name = "Blues"))(length(breaksList))
  }
  pheatmap::pheatmap(m1, color = colours, cutree_rows = 3, cutree_cols = 3)
  
  if(print){
    return(present_terms)
  }
  
}

printGO <- function(x){
  for(go in x){
    print(GOTERM[[go]])
    cat("--------------------------------------\n")
  }
}

# gophr_broke
#'
#' Plot GO information for a list of genes
#' @param genes A list of gene names [geneSymbol]
#' @param burrow Print GO descriptions for each gene in list [binary]
#' @param tunnel Print the list of GO terms [binary]
#' @param waffle Limit the number of terms/desriptions to top N [int]
go_home <- function(genes, burrow=FALSE, waffle=NULL, tunnel=FALSE, go_class=c("GO:BP", "GO:MF", "GO:CC")){
  genes <- unique(genes)
  genes <- matrix(genes)
  
  ontologies <- goseq::getgo(genes, "dm3", "geneSymbol", fetch.cats=go_class)
  
  for(i in 1:length(ontologies)){
    if(is.na(names(ontologies[i])))
      next
    if(missing(waffle)){
      go_terms <- ontologies[[i]]
    }
    else{
      go_terms <- ontologies[[i]][1:waffle]
    }
    if(tunnel){
      cat("Gene: ", names(ontologies[i]), "\n")
      cat("Mappings: ", go_terms, "\n") 
    }
    if(burrow){
      cat("Gene: ", names(ontologies[i]), "\n")
      printGO(go_terms)
    }
  }
  
}





get_down <- function(go_terms, go_class='BP', direction='over', slim_file = 'goslim_generic.obo', evidenceCode="TAS"){
  myCollection <- GOCollection(go_terms)
  
  slim <- getOBOCollection(slim_file)
  slim_counts <- goSlim(myCollection, slim, go_class, evidenceCode %in% evidenceCode)
  slim_counts <- slim_counts %>% 
    dplyr::filter(Count > 0)
  
  
  fill_col <- '#3B8FC7'
  if(direction=='under'){
    fill_col <- '#BF3131'
  }
  
  title <- paste0(direction, "-represented ", go_class, " terms")
  
  p <- ggplot(slim_counts)
  p <- p + geom_bar(aes(fct_reorder(Term, Percent), Percent, fill=fill_col), stat='identity')
  p <- p + scale_x_discrete("Percentage")
  p <- p + coord_flip()
  p <- p + theme(axis.title.y = element_blank()) + cleanTheme()
  p <- p + scale_fill_identity()
  p <- p + ggtitle(title)
  p
}











# hm_1 <- function(combined=each_sample){
#   combined$term <- substr(combined$term, start = 1, stop = 50)
#   df <- melt(combined, id='term')
#   df_f <- df %>% 
#     dplyr::group_by(term) %>% 
#     dplyr::filter(!is.na(value)) %>% 
#     dplyr::add_tally()
#   
#   df2 <- ddply(df_f, .(variable), transform,
#                rescale = rescale(-value))
#   
#   df2 <- df2 %>% 
#     # dplyr::filter(rescale > 0.1) %>% 
#     dplyr::filter(n > 1) %>% 
#     droplevels()
#   
#   base_size <- 9
#   
#   textcol <- "grey40"
#   
#   p <- ggplot(df2, aes(variable, term))
#   p <- p + geom_tile(aes(fill = rescale), colour="white",size=0.25, show.legend=FALSE)
#   
#   p <- p + scale_fill_gradient(low = "#36648B", high = "#87CEFF", na.value="grey90")
#   # p <- p + scale_fill_manual(values=c("#d53e4f","#f46d43","#fdae61",
#   # "#fee08b","#e6f598","#abdda4","#ddf1da"), na.value="grey90")
#   p <- p + scale_x_discrete(expand = c(0, 0))
#   p <- p + scale_y_discrete(expand = c(0, 0))
#   p <- p + theme_grey(base_size = base_size)
#   p <- p + theme(axis.text.x = element_text(size=10,colour=textcol,
#                                             angle = 90, hjust = 1, vjust = 0.5),
#                  axis.title.x = element_blank(),
#                  axis.title.y = element_blank(),
#                  axis.text.y = element_text(size=15,vjust = 0.2,colour=textcol)
#   )
#   p <- p + cleanTheme()
#   p
# }
# hm_2 <- function(combined=all_samples, out_file='all_samples_heatmap', write=NULL, p_threshold = 0.01, go_class=c("BP", "MF", "CC")){
#   
#   scaled <- combined %>% 
#     dplyr::filter(ontology %in% go_class) %>% 
#     dplyr::mutate(term = substr(term, start = 1, stop = 50)) %>% 
#     dplyr::mutate_at(rescale, .vars = vars(-term, -ontology)) %>% 
#     droplevels()
#   
#   df <- melt(combined, id='term') %>% 
#     df_f <- df %>% 
#     dplyr::group_by(term) %>% 
#     dplyr::filter(!is.na(value)) %>% 
#     dplyr::add_tally()
#   
#   
#   scale_this <- function(x) as.vector(scale(x))
#   new <- df_f %>% 
#     dplyr::mutate(new = scale(-value))
#   
#   
#   df2 <- plyr::ddply(df_f, .(variable), transform,
#                      rescale = rescale(-value))
#   
#   df2 <- df2 %>% 
#     # dplyr::filter(rescale > 0.1) %>% 
#     dplyr::filter(n > 1) %>% 
#     dplyr::select(-n, -value) %>% 
#     droplevels()
#   
#   
#   # df3 <- spread(df2, key = variable, value = rescale)
#   
#   # df3 <- dcast(df2, term ~ variable)
#   df4 <- reshape(df2, idvar = "term", timevar = "variable", direction = "wide")
#   
#   colnames(df4) <- gsub(x = names(df4), pattern = "rescale.", replacement = "")  
#   
#   # colnames(df4) <- c("term", levels(df2$variable))
#   
#   df5 <- df4
#   df5$term <- NULL
#   m1 <- as.matrix(df5)
#   rownames(m1) <- df4$term
#   m1[is.na(m1)] <- 0
#   
#   # heatmap(m1, Colv = NA, Rowv = NA, scale="column")
#   
#   if(!missing(write)){
#     png(paste0(out_file, '.png'),    # create PNG for the heat map        
#         width = 6*400,        # 5 x 300 pixels
#         height = 4*400,
#         res = 300,            # 300 pixels per inch
#         pointsize = 7)        # smaller font size
#   }
#   
#   heatmap(m1, scale="column", cexRow=1.5, col=colorRampPalette(brewer.pal(8, "Blues"))(25))
#   dev.off()
#   heatmap(m1, scale="column", cexRow=1.5, col=colorRampPalette(brewer.pal(8, "Blues"))(25))
#   
# }
# 
# 
# # add variable: performance
# set.seed(41)
# df.team_data$performance <- rnorm(nrow(df.team_data))
# 
# #inspect
# head(df.team_data)
# 
# ggplot(data = df.team_data, aes(x = metrics, y = teams)) +
#   geom_tile(aes(fill = performance)) 
# 
# # wide <- reshape(long, idvar = "sample", timevar = "term", direction = "wide")
# 
# df <- data_frame(term = c("term1", "term2", "term3"), sample = c("S1", "S2", "S3"))
# 
# go_matrix <- as.matrix(GO_wall$over_represented_pvalue)
# row.names(go_matrix) <- GO_wall$term
# colnames(go_matrix) <- c('sample1', 'sample2')
# go_matrix$sample2 <- go_matrix$sample1
# View(go_matrix)
# heatmap(go_matrix)
# 
# adjusted_hits=GO_wall$category[p.adjust(GO_wall$over_represented_pvalue,method="BH")<.05]
# View(adjusted_hits)
# 
# library(GO.db)
# for(go in adjusted_hits){
#   print(GOTERM[[go]])
#   cat("--------------------------------------\n")
# }
# 
# KEGG=goseq(pwf,"dm3","geneSymbol",test.cats="KEGG")
# View(KEGG)
# 
# 
# go_map <- function(mat){
#   library(scales)
#   nba <- read.csv("http://datasets.flowingdata.com/ppg2008.csv")
#   nba$Name <- with(nba, reorder(Name, PTS))
#   
#   nba.m <- melt(nba)
#   
#   
#   nba.m <- ddply(nba.m, .(variable), transform,
#                  rescale = rescale(value))
#   
#   p <- ggplot(nba.m, aes(variable, Name))
#   p <- p + geom_tile(aes(fill = rescale), colour = "white")
#   p <- p + scale_fill_gradient(low = "white", high = "steelblue")
#   p
# }