#### Transcriptomics analysis - TMBroiler Muscle Project
#### Shadi Shahatit - RA, JUST 2024-2025
# libraries ---------------------------------------------------------------

sys_dir <- "/home/shadi/Desktop/transcriptome_alzghoul_SS/"

library(ggrepel)
library(DEGreport)
library(ggiraph)
library(ggtree)
library(gdtools)
library('yulab.utils')
library(clusterProfiler)
library(topGO)
library(GO.db)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridisLite)
library(simplifyEnrichment)
library(biomaRt)
library(org.Gg.eg.db)
library(tidyverse)    # dplyr, tidyr, ggplot2, readr
library(forcats)
library(stringr)
library(ggnewscale)
library(readxl)
library(writexl)
library(tximport)
library(txdbmaker)
library(AnnotationDbi)
library(clusterProfiler)
library(DESeq2)
library(ggsci)
library(patchwork)
library(gridExtra)
library(UpSetR)
library(ggVennDiagram)
library(ggvenn)
library(VennDiagram)
library(viridis)
library(enrichplot)
library(tibble)

# Gene-level quantification - tximport -----------------------------------

## tximport

## tximport to get gene-level quantification from Salmon transcript-level output

## convert GTF to TxDb

gtf_dir <- file.path(sys_dir, "TMBroilers_Transcriptomics/Ggallus_Refs/gtf_r115/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.115.chr.gtf.gz")
txdb_Ggallus <- txdbmaker::makeTxDbFromGFF(gtf_dir,
                                           format="gtf",
                                           dataSource="Ensembl",
                                           organism="Gallus_gallus",
                                           taxonomyId=9031)

## extract transcript-to-gene mapping from TxDb object

k <- keys(txdb_Ggallus, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb_Ggallus, k, "GENEID", "TXNAME")

length(unique(tx2gene$GENEID))
length(unique(tx2gene$TXNAME))
keys(txdb_Ggallus, keytype = "GENEID")

## upload salmon for tximport and sum into gene level quants

quants_dir <- file.path(sys_dir, "TMBroilers_Transcriptomics/muscle/quant_r115")
files <- list.files(quants_dir, pattern = "quant.sf", recursive = TRUE, full.names = TRUE)
names(files) <- basename(dirname(files))
TxiSalmon_muscle <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = T)
# txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = T)
# txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = T)

head(TxiSalmon_muscle$counts)     # raw_counts
head(TxiSalmon_muscle$abundance)  # TPM
head(TxiSalmon_muscle$length)     # gene_lengths (bps)
ncol(TxiSalmon_muscle$counts)     # 25 is the number of our samples for muscle 

## get exp matrix at the level of transcripts as well

TxiTranscripts_muscle <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE,
                                  txOut = TRUE)
transcript_counts_muscle <- as.data.frame(TxiTranscripts_muscle$counts)

# to remove transcript versions
transcript_counts_muscle$TXNAME <- gsub("\\.\\d+$", "", rownames(transcript_counts_muscle))
transcript_counts_muscle_mapped <- merge(transcript_counts_muscle, tx2gene, by = "TXNAME")
transcripts_per_gene_muscle <- transcript_counts_muscle_mapped %>%
  dplyr::select(TXNAME, GENEID) %>%
  distinct() %>%
  group_by(GENEID) %>%
  summarise(transcript_count = n())

## plot transcript counts per gene

ggplot(subset(transcripts_per_gene_muscle,transcript_count<100), aes(x = transcript_count)) +
  geom_histogram(binwidth = 1, fill = "#323232") +
  labs(title = "", x = "Transcripts per gene", y = "Gene counts") +
  theme_classic()

summary(transcripts_per_gene_muscle$transcript_count)
mean(transcripts_per_gene_muscle$transcript_count)
sd(transcripts_per_gene_muscle$transcript_count)
t.test(transcripts_per_gene_muscle$transcript_count, conf.level = 0.95)$conf.int
median(transcripts_per_gene_muscle$transcript_count)
IQR(transcripts_per_gene_muscle$transcript_count)
quantile(transcripts_per_gene_muscle$transcript_count, probs = c(0.25, 0.75))
nrow(transcripts_per_gene_muscle)
sum(transcripts_per_gene_muscle$transcript_count)

# Metadata and ExpMat preparation  ---------------------------------------------------

## prep your metadata - muscle

## create a table with consistent naming with txi salmon files

colnames(TxiSalmon_muscle$counts)

metadata_muscle <- tibble(sample_id = colnames(TxiSalmon_muscle$counts)) %>%
  mutate(id = sample_id,
         # id = str_remove(sample_id, "_quant$"),
         lower = tolower(sample_id),
         tissue = "muscle",
         day = case_when(
           str_detect(lower, "d22|a22|at22|ct22") ~ "D22",
           str_detect(lower, "d19") ~ "D19",
           str_detect(lower, "(^|_)p7|(^|_)d7") ~ "D7",
           TRUE ~ NA_character_),
         treatment_lvl1 = case_when(
           str_detect(lower, "(^|_)con(_|$)") ~ "none",
           str_detect(lower, "(^|_)tm(_|$)") ~ "TM",
           TRUE ~ NA_character_),
         treatment_lvl2 = case_when(
           str_detect(lower, "^a") ~ "acute_HS",
           str_detect(lower, "^c") ~ "chronic_HS",
           TRUE ~ "none")) %>%
  transmute(id,
            sample_id,
            tissue,
            day = factor(day, levels = c("D19","D7","D22")),
            treatment_lvl1 = factor(treatment_lvl1, levels = c("none","TM")),
            treatment_lvl2 = factor(treatment_lvl2, levels = c("none","acute_HS","chronic_HS")),
            condition = factor(case_when(
              (treatment_lvl1 %in% c("none", NA)) & (treatment_lvl2 %in% c("none", NA)) ~ "ctrl",
              (treatment_lvl1 %in% c("none", NA)) & !(treatment_lvl2 %in% c("none", NA)) ~ as.character(treatment_lvl2),
              !(treatment_lvl1 %in% c("none", NA)) & (treatment_lvl2 %in% c("none", NA)) ~ as.character(treatment_lvl1),
              TRUE ~ paste(treatment_lvl1, treatment_lvl2, sep = "_")
            ), levels = c("ctrl","TM","acute_HS","chronic_HS","TM_acute_HS","TM_chronic_HS"))) %>%
  column_to_rownames("sample_id")

stopifnot(identical(rownames(metadata_muscle), colnames(TxiSalmon_muscle$counts)))

## change the names of conditions

rename_map <- c(
  "ctrl"          = "Con_N",
  "TM"            = "TM_N",
  "chronic_HS"    = "Con_CHS",
  "TM_chronic_HS" = "TM_CHS",
  "acute_HS"      = "Con_AHS",
  "TM_acute_HS"   = "TM_AHS")

metadata_muscle$condition <- recode(metadata_muscle$condition, !!!rename_map)
unique(metadata_muscle$condition)

## since you have two lvls in condition (TM + HS), prep a metadata and ExpMat for each

## TM analysis

metadata_muscle_TM <- metadata_muscle %>%
  filter(condition == "Con_N" | condition == "TM_N")
# metadata_muscle_TM <- metadata_muscle %>%
#   filter(!grepl("HS", condition))
keep_TMonly <- metadata_muscle_TM$id
# TxiSalmon_muscle_TM <- TxiSalmon_muscle$counts[, keep_TMonly]
TxiSalmon_muscle_TM <- list(
  counts    = TxiSalmon_muscle$counts[, keep_TMonly, drop = FALSE],
  abundance = TxiSalmon_muscle$abundance[, keep_TMonly, drop = FALSE],
  length    = TxiSalmon_muscle$length[, keep_TMonly, drop = FALSE],
  countsFromAbundance = TxiSalmon_muscle$countsFromAbundance[])
(identical(rownames(metadata_muscle_TM), colnames(TxiSalmon_muscle_TM$counts)))
(identical(rownames(metadata_muscle_TM), colnames(TxiSalmon_muscle$counts)))

## HS analysis

metadata_muscle_D22 <- metadata_muscle %>%
  filter(day == "D22")
keep_D22only <- metadata_muscle_D22$id
# TxiSalmon_muscle_D22 <- TxiSalmon_muscle$counts[, keep_D22only]
TxiSalmon_muscle_D22 <- list(
  counts    = TxiSalmon_muscle$counts[, keep_D22only, drop = FALSE],
  abundance = TxiSalmon_muscle$abundance[, keep_D22only, drop = FALSE],
  length    = TxiSalmon_muscle$length[, keep_D22only, drop = FALSE],
  countsFromAbundance = TxiSalmon_muscle$countsFromAbundance[])
(identical(rownames(metadata_muscle_D22), colnames(TxiSalmon_muscle_D22$counts)))
(identical(rownames(metadata_muscle_D22), colnames(TxiSalmon_muscle$counts)))

# DEseq2 modeling: HS analysis --------------------------------------------

## (1) Per day/condition pairwise comp      - Wald

## HS analysis

## create a DESeq object
dds_muscle_D22_00 <- DESeqDataSetFromTximport(
  txi = TxiSalmon_muscle_D22,
  colData = metadata_muscle_D22,
  design = ~ condition) # Pairwise comp of all condition at D22; Wald

## Pre-filter low-expression genes
dds_muscle_D22_01 <- dds_muscle_D22_00[rowSums(counts(dds_muscle_D22_00)) >= 10, ]
# dds$day <- droplevels(dds$day)
# dds$condition <- droplevels(dds$condition)

## DESeq modeling
dds_muscle_D22 <- DESeq(dds_muscle_D22_01)

head(counts(dds_muscle_D22))
colData(dds_muscle_D22)
summary(dds_muscle_D22)
levels(dds_muscle_D22$condition)
sizeFactors(dds_muscle_D22)
dispersions(dds_muscle_D22)

# Perform DESeq with LRT for global analysis
# dds_MA <- DESeq(dds_MA, test = "LRT", reduced = ~ 1)

## extract results (AEG_full_info_muscle_D22; DEG_full_info_muscle_D22)

padj_thresh <- 0.05
lfc_thresh  <- 1

## New way: define res and no need to flip lcf

pairs_target <- list(
  c("TM_N","Con_N"),
  c("Con_AHS","Con_N"),
  c("Con_AHS","TM_N"),
  c("TM_AHS","Con_N"),
  c("TM_AHS","TM_N"),
  c("TM_AHS","Con_AHS"),
  c("Con_CHS","Con_N"),
  c("Con_CHS","TM_N"),
  c("TM_CHS","Con_N"),
  c("TM_CHS","TM_N"),
  c("TM_CHS","Con_CHS"),
  c("Con_CHS","Con_AHS"),
  c("TM_AHS","Con_CHS"),
  c("TM_CHS","Con_AHS"),
  c("TM_CHS","TM_AHS"))

AEG_full_info_muscle_D22_list <- lapply(seq_along(pairs_target), function(i){
  g1 <- pairs_target[[i]][1]; g2 <- pairs_target[[i]][2]
  labs <- paste(g1, "vs", g2)
  res <- try(results(dds_muscle_D22, contrast = c("condition", g1, g2)), silent = TRUE)
  tibble::as_tibble(res, rownames = "Gene") %>% mutate(contrast = labs)
})

AEG_full_info_muscle_D22 <- dplyr::bind_rows(AEG_full_info_muscle_D22_list)

AEG_full_info_muscle_D22$contrast <- factor(AEG_full_info_muscle_D22$contrast,
                                            levels = vapply(pairs_target, function(p) paste(p[1],"vs",p[2]), ""))

DEG_full_info_muscle_D22 <- AEG_full_info_muscle_D22 %>%
  filter(!is.na(padj),
         padj < padj_thresh,
         abs(log2FoldChange) > lfc_thresh) %>%
  mutate(direction = if_else(log2FoldChange > 0, "upregulated", "downregulated"))

length(unique(AEG_full_info_muscle_D22$Gene)) < length(unique(rownames(TxiSalmon_muscle_D22$counts)))
nrow(DEG_full_info_muscle_D22)

## append different gene ids and download

collapse_unique <- function(x) {
  x <- unique(na.omit(x))
  if (length(x) == 0) NA_character_ else paste(x, collapse = ";")
}

IDmap_Gg <- AnnotationDbi::select(
  org.Gg.eg.db,
  keys    = unique(DEG_full_info_muscle_D22$Gene),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL","SYMBOL","ENTREZID","REFSEQ")) %>%
  as_tibble() %>%
  group_by(ENSEMBL) %>%
  summarise(
    SYMBOL   = collapse_unique(SYMBOL),
    ENTREZID = collapse_unique(ENTREZID),
    REFSEQ   = collapse_unique(REFSEQ),
    .groups  = "drop")

DEG_full_info_muscle_D22_all_anno <- DEG_full_info_muscle_D22 %>%
  left_join(IDmap_Gg, by = c(Gene = "ENSEMBL")) %>%
  mutate(SYMBOL = if_else(is.na(SYMBOL) | SYMBOL == "", Gene, SYMBOL))

write_xlsx(
  DEG_full_info_muscle_D22_all_anno,
  file.path(sys_dir, "TMBroilers_Transcriptomics", "TMomics_DEG_muscle_D22_allanno_allcomps_v5_final.xlsx"))

# DEGs visualizations (Bar + Volcano + Venn diagrams + Upset plots): HS analysis -----------------------------------------------------------

## Bar plots

## count and pivot DEGs per contrast

DEG_muscle_D22_sum <- DEG_full_info_muscle_D22 %>%
  group_by(contrast, direction) %>%
  summarize(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = direction, values_from = count, values_fill = 0) %>%
  mutate(total = upregulated + downregulated)

DEG_muscle_D22_sum_long <- DEG_muscle_D22_sum %>%
  pivot_longer(cols = c("upregulated", "downregulated", "total"),
               names_to = "category", values_to = "count") %>% 
  mutate(category = factor(category, levels = c("downregulated", "upregulated", "total")))

ggplot(subset(DEG_muscle_D22_sum_long), aes(x = contrast, y = count, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = count), position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5) +
  theme_classic() +
  labs(title = "",
       x = "Conditions",
       y = "DEGs counts",
       fill = "DEGs category") +
  scale_fill_manual(values = c("#FFD662FF", "#00539CFF", "#949398FF")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
# theme(
#   plot.title = element_text(size = 16, face = "bold"),
#   axis.title.x = element_text(size = 14, face = "bold"),
#   axis.title.y = element_text(size = 14, face = "bold"))

## Volcano plots - v2

collapse_unique <- function(x) {
  x <- unique(na.omit(x))
  if (length(x) == 0) NA_character_ else paste(x, collapse = ";")
}

IDmap_Gg <- AnnotationDbi::select(
  org.Gg.eg.db,
  keys    = unique(DEG_full_info_muscle_D22$Gene),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL","SYMBOL")) %>%
  as_tibble() %>%
  group_by(ENSEMBL) %>%
  summarise(
    SYMBOL   = collapse_unique(SYMBOL),
    .groups  = "drop")

DEG_full_info_muscle_D22_all_anno <- DEG_full_info_muscle_D22 %>%
  left_join(IDmap_Gg, by = c(Gene = "ENSEMBL")) %>%
  # group_by(contrast,direction) %>%
  # slice_max(order_by = log2FoldChange, n = 10) %>%
  mutate(SYMBOL = if_else(is.na(SYMBOL) | SYMBOL == "" | grepl(";", SYMBOL), NA, SYMBOL))

volc_cols <- c("Upregulated" = "#22A884FF", "Downregulated" = "#440154FF", "Not Significant" = "#B3B3B3")

VolcanoPlot_label <- function(df, df_DEG, ttl){
  df2 <- df %>% mutate(direction = case_when(
    log2FoldChange >  lfc_thresh & padj < padj_thresh ~ "Upregulated",
    log2FoldChange < -lfc_thresh & padj < padj_thresh ~ "Downregulated",
    TRUE ~ "Not Significant")) %>%
    left_join(df_DEG %>% select(Gene, contrast, SYMBOL), by = c("Gene", "contrast"))
  
  # df_labels <- df2 %>% filter(direction != "Not Significant")
  df_labels <- df2 %>% 
    filter(direction != "Not Significant") %>%
    group_by(contrast, direction) %>%
    slice_max(order_by = log2FoldChange, n = 10) %>%
    ungroup()
  
  upc <- sum(df2$direction=="Upregulated", na.rm=TRUE)
  dnc <- sum(df2$direction=="Downregulated", na.rm=TRUE)
  ymax <- max(-log10(df2$padj), na.rm=TRUE); if(!is.finite(ymax)) ymax <- 1
  
  ggplot(df2, aes(log2FoldChange, -log10(padj), color=direction)) +
    geom_point(size=1.6, alpha=0.8) +
    scale_color_manual(values=volc_cols) +
    geom_hline(yintercept=-log10(padj_thresh), linetype="dashed") +
    geom_vline(xintercept=c(-lfc_thresh, lfc_thresh), linetype="dashed") +
    annotate("text", x=-5.5, y=ymax, label=dnc, size=3.8, fontface="bold") +
    annotate("text", x= 5.5, y=ymax, label=upc, size=3.8, fontface="bold") +
    xlim(-6,6) + labs(title=ttl, x="log2 FC", y="-log10(q)") +
    geom_text_repel(data = df_labels, aes(label = SYMBOL, color = direction),
                    size = 3, max.overlaps = Inf,
                    box.padding = 0.5, point.padding = 0.25,
                    show.legend = FALSE) +
    theme_minimal(base_size=11) +
    theme(panel.grid=element_blank(),
          panel.border=element_rect(color="black", fill=NA, linewidth=0.6),
          legend.position="none",
          plot.title=element_text(hjust=0.5, face="bold"))
}

AEG_full_info_muscle_D22_set1 <- AEG_full_info_muscle_D22 %>%
  filter(contrast %in% GeneSet_PairSet_target$set1) %>%
  droplevels()

AEG_full_info_muscle_D22_set2 <- AEG_full_info_muscle_D22 %>%
  filter(contrast %in% GeneSet_PairSet_target$set2) %>%
  droplevels()

AEG_full_info_muscle_D22_set3 <- AEG_full_info_muscle_D22 %>%
  filter(contrast %in% "Con_CHS vs Con_AHS") %>%
  droplevels()

## set 1

VolcPlot_muscle_D22_set1 <- lapply(
  split(AEG_full_info_muscle_D22_set1, AEG_full_info_muscle_D22_set1$contrast),
  function(d) VolcanoPlot_label(d %>% select(Gene, log2FoldChange, padj, contrast),
                                DEG_full_info_muscle_D22_all_anno,
                                ttl = unique(d$contrast)[1]))

wrap_plots(VolcPlot_muscle_D22_set1, ncol = 3) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20, face = "bold", family = "Times New Roman",color = "white"))

## set 2

VolcPlot_muscle_D22_set2 <- lapply(
  split(AEG_full_info_muscle_D22_set2, AEG_full_info_muscle_D22_set2$contrast),
  function(d) VolcanoPlot_label(d %>% select(Gene, log2FoldChange, padj, contrast),
                                DEG_full_info_muscle_D22_all_anno,
                                ttl = unique(d$contrast)[1]))

wrap_plots(VolcPlot_muscle_D22_set2, ncol = 3) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20, face = "bold", family = "Times New Roman",color = "white"))

## set 3

VolcPlot_muscle_D22_set3 <- lapply(
  split(AEG_full_info_muscle_D22_set3, AEG_full_info_muscle_D22_set3$contrast),
  function(d) VolcanoPlot_label(d %>% select(Gene, log2FoldChange, padj, contrast),
                                DEG_full_info_muscle_D22_all_anno,
                                ttl = unique(d$contrast)[1]))

wrap_plots(VolcPlot_muscle_D22_set3, ncol = 1) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20, face = "bold", family = "Times New Roman",color = "white"))

## Upset plots + Venn diagrams

GeneSet_muscle_D22 <- DEG_full_info_muscle_D22_all_anno %>%
  group_by(contrast) %>%
  summarize(Gene_Symbols = list(unique(Gene)), .groups = "drop") %>%
  deframe()

GeneSet_PairSet_target <- list(
  set1 = c(
    "TM_N vs Con_N",
    "Con_AHS vs Con_N",
    # "Con_AHS vs TM_N",
    "TM_AHS vs Con_N",
    # "TM_AHS vs TM_N",
    "TM_AHS vs Con_AHS"
  ),
  set2 = c(
    "TM_N vs Con_N",
    "Con_CHS vs Con_N",
    # "Con_CHS vs TM_N",
    "TM_CHS vs Con_N",
    # "TM_CHS vs TM_N",
    "TM_CHS vs Con_CHS"
  ),
  set3 = c(
    "Con_CHS vs Con_AHS",
    "TM_AHS vs Con_CHS",
    "TM_CHS vs Con_AHS",
    "TM_CHS vs TM_AHS"
  ))

GeneSet_muscle_D22_1 <- GeneSet_muscle_D22[GeneSet_PairSet_target$set1]
GeneSet_muscle_D22_2 <- GeneSet_muscle_D22[GeneSet_PairSet_target$set2]

## Upset

GeneSet_Matrix_muscle_D22_1 <- fromList(GeneSet_muscle_D22_1)
GeneSet_Matrix_muscle_D22_2 <- fromList(GeneSet_muscle_D22_2)

UpSetR::upset(GeneSet_Matrix_muscle_D22_1,
              main.bar.color = "#404080",
              sets.bar.color = "#69b3a2",
              order.by = "freq",
              sets = names(GeneSet_muscle_D22_1),
              set_size.show = TRUE)

UpSetR::upset(GeneSet_Matrix_muscle_D22_2,
              main.bar.color = "#404080",
              sets.bar.color = "#69b3a2",
              order.by = "freq",
              sets = names(GeneSet_muscle_D22_2),
              set_size.show = TRUE)

## Venn diagrams

ggVennDiagram(GeneSet_muscle_D22_1, label = "none",
              set_color = viridis(length(GeneSet_muscle_D22_1)))
ggVennDiagram(GeneSet_muscle_D22_2, label = "none",
              set_color = viridis(length(GeneSet_muscle_D22_2)))

ggvenn(
  GeneSet_muscle_D22_1, 
  fill_color = viridis(length(names(GeneSet_muscle_D22_1))),
  # fill_color = c("#440154FF","#22A884FF","#2A788EFF","#FDE725FF"),
  # fill_color = c("#440154FF","#22A884FF","#2A788EFF","#FDE725FF"),
  show_percentage = FALSE,
  stroke_size = 0.5, 
  set_name_size = 4, 
  text_size = 4)

ggvenn(
  GeneSet_muscle_D22_2, 
  fill_color = viridis(length(names(GeneSet_muscle_D22_2))),
  # fill_color = c("#440154FF","#22A884FF","#2A788EFF","#FDE725FF"),
  # fill_color = c("#440154FF","#22A884FF","#2A788EFF","#FDE725FF"),
  show_percentage = FALSE,
  stroke_size = 0.5, 
  set_name_size = 4, 
  text_size = 4)

# Functional enrichment v2 - HS analysis ---------------------------------------------------------------------

## using these dfs
## AEG_full_info_muscle_D22
## DEG_full_info_muscle_D22
## dds_muscle_D22

## using helper functions to perform analysis and plot
## perform_ORA_GO_BP
## perform_ORA_KEGG
## perform_GSEA_GO
## perform_GSEA_KEGG
## ORA_PublishFEBarplot
## GSEA_PublishDotplot

## over-representation analysis ORA

## helper function for ordered panel FC bar plots: ORA table (Description, Fold Enrichment, p-value, contrast)
ORA_PublishFEBarplot <- function(ora_df, n_top = 10, ncol = 3, title_col = "contrast", order_levels = NULL) {
  dfp <- ora_df %>%
    mutate(
      log10p = -log10(p.adjust),
      Desc   = str_wrap(Description, width = 40),
      !!title_col := if (is.null(order_levels)) .data[[title_col]]
      else factor(.data[[title_col]], levels = order_levels)
    ) %>%
    group_by(.data[[title_col]]) %>%
    slice_min(order_by = p.adjust, n = n_top, with_ties = FALSE) %>%
    ungroup()
  
  levs <- if (is.null(order_levels)) levels(factor(dfp[[title_col]])) else order_levels
  plots <- lapply(levs, function(lbl) {
    dd <- dfp %>% filter(.data[[title_col]] == lbl)
    ggplot(dd, aes(x = FoldEnrichment, y = fct_reorder(Desc, FoldEnrichment), fill = log10p)) +
      geom_bar(stat = "identity", width = 0.7) +
      # scale_fill_gradient(low = "#BCE388FF", high = "#22A884FF", name = expression(-log[10](q))) +
      scale_fill_gradient(low = "firebrick", high = "darkgreen", name = expression(-log[10](q))) +
      labs(x = "Fold enrichment", y = NULL, title = lbl) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.y        = element_text(size = 11),
        axis.text.x        = element_text(size = 11),
        axis.title.x       = element_text(size = 13),
        plot.title         = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title       = element_text(size = 13),
        legend.text        = element_text(size = 12),
        panel.grid.major.y = element_blank()
      )
  })
  wrap_plots(plots, ncol = ncol)
}

## DEGs list per contrast with ENSEMBL IDs and a specific pairwise sets

degs_sets_muscle_D22 <-
  DEG_full_info_muscle_D22 %>%
  distinct(contrast, Gene) %>%
  group_by(contrast) %>%
  summarise(genes = list(unique(Gene)), .groups = "drop") %>%
  { setNames(.$genes, .$contrast) }

## in case of targeted pairwise comp

GeneSet_PairSet_target <- list(
  set1 = c(
    "TM_N vs Con_N",
    "Con_AHS vs Con_N",
    "Con_AHS vs TM_N",
    "TM_AHS vs Con_N",
    "TM_AHS vs TM_N",
    "TM_AHS vs Con_AHS"
  ),
  set2 = c(
    # "TM_N vs Con_N",
    "Con_CHS vs Con_N",
    "Con_CHS vs TM_N",
    "TM_CHS vs Con_N",
    "TM_CHS vs TM_N",
    "TM_CHS vs Con_CHS"
  ),
  set3 = c(
    "Con_CHS vs Con_AHS",
    "TM_AHS vs Con_CHS",
    "TM_CHS vs Con_AHS",
    "TM_CHS vs TM_AHS"
  ))

## Background universe genes with ENSEMBL IDs

background_genes_muscle_D22 <- rownames(dds_muscle_D22)

## over-representation analysis ORA (GO:BP)

perform_ORA_GO_BP <- function(gene_set, label) {
  enrichGO(
    gene          = gene_set,
    OrgDb         = org.Gg.eg.db,
    keyType       = "ENSEMBL",
    ont           = "BP",
    universe      = background_genes_muscle_D22,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = FALSE
  ) |>
    as.data.frame() |>
    dplyr::mutate(contrast = label)
}

ORA_GO_muscle_D22_list <- lapply(names(degs_sets_muscle_D22), function(lbl) {
  perform_ORA_GO_BP(degs_sets_muscle_D22[[lbl]], lbl)
})

ORA_GO_muscle_D22 <- bind_rows(ORA_GO_muscle_D22_list) %>%
  mutate(level = "all_DEGs") 

ORA_GO_muscle_D22_PlotmePanel <- ORA_PublishFEBarplot(ORA_GO_muscle_D22, n_top = 10, ncol = 3, title_col = "contrast")

ORA_PublishFEBarplot(ORA_GO_muscle_D22, n_top = 10, ncol = 3, title_col = "contrast",
                     order_levels = GeneSet_PairSet_target$set1)
ORA_PublishFEBarplot(ORA_GO_muscle_D22, n_top = 10, ncol = 3, title_col = "contrast",
                     order_levels = GeneSet_PairSet_target$set2)
ORA_PublishFEBarplot(ORA_GO_muscle_D22, n_top = 10, ncol = 3, title_col = "contrast",
                     order_levels = GeneSet_PairSet_target$set3)

## EXTRA: one can remove redundant enriched GO terms using similarity matrix analysis via simplifyEnrichment::simplify()

sem_bp_gga <- GOSemSim::godata(annoDb = "org.Gg.eg.db", ont = "BP", computeIC = FALSE)

clustermap  <- simplifyEnrichment::simplifyGO(unique(ORA_GO_muscle_D22$ID),
                                              ont = "BP", sem_data = sem_bp_gga, cutoff = 0.7, plot = FALSE)

ORA_GO_muscle_D22_simp <- ORA_GO_muscle_D22 %>%
  dplyr::left_join(clustermap, by = c("ID" = "id")) %>%
  dplyr::group_by(contrast, cluster) %>%
  dplyr::slice_min(order_by = p.adjust, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(level = "all_DEGs")

ORA_GO_muscle_D22_PlotmePanel_simp <- ORA_PublishFEBarplot(ORA_GO_muscle_D22_simp, n_top = 10, ncol = 3, title_col = "contrast")

ORA_PublishFEBarplot(ORA_GO_muscle_D22_simp, n_top = 10, ncol = 3, title_col = "contrast",
                     order_levels = GeneSet_PairSet_target$set1)

## gene set enrichment analysis GSEA

## v2

GeneSet_PairSet_target <- list(
  set1 = c(
    "TM_N vs Con_N",
    "Con_AHS vs Con_N",
    "Con_AHS vs TM_N",
    "TM_AHS vs Con_N",
    "TM_AHS vs TM_N",
    "TM_AHS vs Con_AHS"
  ),
  set2 = c(
    "TM_N vs Con_N",
    "Con_CHS vs Con_N",
    "Con_CHS vs TM_N",
    "TM_CHS vs Con_N",
    "TM_CHS vs TM_N",
    "TM_CHS vs Con_CHS"
  ),
  set3 = c(
    "Con_CHS vs Con_AHS",
    "TM_AHS vs Con_CHS",
    "TM_CHS vs Con_AHS",
    "TM_CHS vs TM_AHS"
  ))

## helper function for GSEA dot plots: ORA table (Description, NES, p.adjust, contrast; optional ONTOLOGY)
GSEA_PublishDotplot <- function(gsea_df, n_top = 10, ncol = 3, order_levels = NULL) {
  
  ## drop rows with NA NES and p.adjust
  gsea_df <- gsea_df %>%
    dplyr::filter(!is.na(NES), !is.na(p.adjust))
  
  df <- gsea_df %>%
    dplyr::mutate(
      log10padj = -log10(p.adjust),
      Desc      = stringr::str_wrap(Description, 40),
      contrast  = if (is.null(order_levels)) contrast else factor(contrast, levels = order_levels),
      ## better define sign and avoid NA facet
      sign      = dplyr::case_when(
        NES > 0  ~ "up",
        NES < 0  ~ "down",
        TRUE     ~ NA_character_
      )
    ) %>%
    dplyr::group_by(contrast, sign) %>%
    dplyr::slice_min(order_by = p.adjust, n = n_top, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  if (!nrow(df)) return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle("No enriched sets"))
  
  ## only plot contrasts that actually exist after filtering
  levs  <- unique(as.character(df$contrast))
  
  plots <- lapply(levs, function(lbl){
    dd <- df[df$contrast == lbl, , drop = FALSE]
    
    max_abs <- max(abs(dd$NES), na.rm = TRUE)
    # ggplot2::ggplot(dd,
    #                 ggplot2::aes(x = NES,
    #                              y = forcats::fct_reorder(Desc, NES),
    #                              size = log10padj,
    #                              color = sign)) +
    #   ggplot2::geom_point(alpha = 0.9) +
    #   ggplot2::scale_color_manual(values = c(down = "#3B4CC0", up = "#B40426")) +
    #   ggplot2::scale_x_continuous(limits = c(-max_abs-1, max_abs+1),
    #                               expand = ggplot2::expansion(mult = c(0.1, 0.1))) +
    #   ggplot2::labs(x = "NES", y = NULL, title = lbl,
    #                 size = expression(-log[10](adj.P))) +
    #   ## theme closer to your example
    #   ggplot2::theme_bw(base_size = 13) +
    #   ggplot2::theme(
    #     panel.grid.major.y = ggplot2::element_blank(),
    #     panel.grid.minor   = ggplot2::element_blank(),
    #     panel.border       = ggplot2::element_rect(color = "black", linewidth = 0.6),
    #     axis.line          = ggplot2::element_line(color = "black"),
    #     plot.title         = ggplot2::element_text(hjust = 0.5, face = "bold"),
    #   ) +
    #   ggplot2::facet_grid(. ~ sign)
    ggplot2::ggplot(dd, ggplot2::aes(x = NES, y = forcats::fct_reorder(Desc, NES), size = log10padj, color = sign)) +
      ggplot2::geom_point(alpha = 0.9) +
      # ggplot2::scale_color_manual(values = c(down = "firebrick", up = "darkgreen")) +
      ggplot2::scale_color_manual(values = c(down = "#3B4CC0", up = "#B40426")) +
      ggplot2::scale_x_continuous(limits = c(-max_abs-1, max_abs+1),
                                  expand = ggplot2::expansion(mult = c(0.1, 0.1))) +
      ggplot2::labs(x = "NES", y = NULL, title = lbl, size = expression(-log[10](adj.P))) +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(
        axis.text.y  = ggplot2::element_text(size = 10),
        axis.text.x  = ggplot2::element_text(size = 10),
        plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.title = ggplot2::element_text(size = 11),
        legend.text  = ggplot2::element_text(size = 10)
      ) +
      ggplot2::facet_grid(. ~ sign)
    # ggplot2::ggplot(dd,
    #                 ggplot2::aes(x = NES,
    #                              y = forcats::fct_reorder(Desc, NES),
    #                              fill = sign)) +
    #   ggplot2::geom_col(width = 0.7) +
    #   ggplot2::scale_fill_manual(values = c(down = "firebrick", up = "darkgreen")) +
    #   ggplot2::scale_fill_manual(values = c(down = "#3B4CC0", up = "#B40426")) +
    #   ggplot2::labs(x = "NES", y = NULL, title = lbl, fill = "direction") +
    #   ggplot2::theme_minimal(base_size = 13) +
    #   ggplot2::theme(
    #     axis.text.y  = ggplot2::element_text(size = 10),
    #     axis.text.x  = ggplot2::element_text(size = 10),
    #     plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold"),
    #     legend.title = ggplot2::element_text(size = 11),
    #     legend.text  = ggplot2::element_text(size = 10)
    #   ) +
    #   ggplot2::facet_grid(. ~ sign)
  })
  
  patchwork::wrap_plots(plots, ncol = ncol)
}

## gene set enrichment analysis GSEA (GO)

## get the ranked gene lists based on expression for each contrast

ranked_aegs_sets_muscle_D22 <- {
  df <- AEG_full_info_muscle_D22[, c("Gene","log2FoldChange","contrast")]
  cn <- unique(df$contrast)
  setNames(lapply(cn, function(x){
    v <- df[df$contrast == x, c("Gene","log2FoldChange")]
    v <- v[is.finite(v$log2FoldChange), , drop = FALSE]
    v <- v[!is.na(v$Gene) & !duplicated(v$Gene), , drop = FALSE]
    g <- v$log2FoldChange; names(g) <- v$Gene
    sort(g, decreasing = TRUE)
  }), cn)
}

## do the GSEA analysis via GO:BP

perform_GSEA_GO <- function(ranked_gene_list, ontology) {
  gseGO(
    geneList      = ranked_gene_list,
    ont           = ontology,
    keyType       = "ENSEMBL",
    minGSSize     = 3,
    maxGSSize     = 800,
    pvalueCutoff  = 0.05,
    verbose       = FALSE,
    OrgDb         = org.Gg.eg.db,
    pAdjustMethod = "BH"
  )
}

ontologies <- c("BP","MF","CC")
GSEA_GO_muscle_D22_list <- setNames(lapply(names(ranked_aegs_sets_muscle_D22), function(cn){
  gl <- ranked_aegs_sets_muscle_D22[[cn]]
  setNames(lapply(ontologies, function(ont) perform_GSEA_GO(gl, ont)), ontologies)
}), names(ranked_aegs_sets_muscle_D22))

GSEA_GO_muscle_D22 <- bind_rows(
  lapply(names(GSEA_GO_muscle_D22_list), function(cn){
    bind_rows(lapply(ontologies, function(ont){
      as.data.frame(GSEA_GO_muscle_D22_list[[cn]][[ont]]) %>%
        dplyr::mutate(contrast = cn, ONTOLOGY = ont)
    }))
  })) %>%
  dplyr::mutate(level = "all_genes")

GSEA_GO_BP_muscle_D22_PlotmePanel  <- GSEA_PublishDotplot(GSEA_GO_muscle_D22[GSEA_GO_muscle_D22$ONTOLOGY == "BP", ],
                                                          n_top = 5, ncol = 3)
GSEA_GO_MF_muscle_D22_PlotmePanel  <- GSEA_PublishDotplot(GSEA_GO_muscle_D22[GSEA_GO_muscle_D22$ONTOLOGY == "MF", ],
                                                          n_top = 10, ncol = 3)
GSEA_GO_CC_muscle_D22_PlotmePanel  <- GSEA_PublishDotplot(GSEA_GO_muscle_D22[GSEA_GO_muscle_D22$ONTOLOGY == "CC", ],
                                                          n_top = 10, ncol = 3)

GSEA_PublishDotplot(GSEA_GO_muscle_D22[GSEA_GO_muscle_D22$ONTOLOGY == "BP", ],
                    n_top = 5, ncol = 3,
                    order_levels = intersect(GeneSet_PairSet_target$set1,
                                             unique(GSEA_GO_muscle_D22$contrast)))

GSEA_GO_muscle_D22_BP_set1a_0.05 <- GSEA_GO_muscle_D22[GSEA_GO_muscle_D22$ONTOLOGY == "BP", ] %>%
  filter(contrast %in% c("TM_N vs Con_N","Con_AHS vs Con_N","Con_AHS vs TM_N","TM_AHS vs Con_N"))
GSEA_GO_muscle_D22_BP_set1b_0.05 <- GSEA_GO_muscle_D22[GSEA_GO_muscle_D22$ONTOLOGY == "BP", ] %>%
  filter(contrast %in% c("TM_AHS vs TM_N","TM_AHS vs Con_AHS"))
GSEA_GO_muscle_D22_BP_set2a_0.05 <- GSEA_GO_muscle_D22[GSEA_GO_muscle_D22$ONTOLOGY == "BP", ] %>%
  filter(contrast %in% c("TM_N vs Con_N","Con_CHS vs Con_N","Con_CHS vs TM_N","TM_CHS vs Con_N"))
GSEA_GO_muscle_D22_BP_set2b_0.05 <- GSEA_GO_muscle_D22[GSEA_GO_muscle_D22$ONTOLOGY == "BP", ] %>%
  filter(contrast %in% c("TM_CHS vs TM_N","TM_CHS vs Con_CHS"))
GSEA_GO_muscle_D22_BP_set3a_0.05 <- GSEA_GO_muscle_D22[GSEA_GO_muscle_D22$ONTOLOGY == "BP", ] %>%
  filter(contrast %in% c("Con_CHS vs Con_AHS"))

GSEA_PublishDotplot(GSEA_GO_muscle_D22_BP_set1a_0.05, n_top = 5, ncol = 2, order_levels = GeneSet_PairSet_target$set1)
GSEA_PublishDotplot(GSEA_GO_muscle_D22_BP_set1b_0.05, n_top = 5, ncol = 2, order_levels = GeneSet_PairSet_target$set1)
GSEA_PublishDotplot(GSEA_GO_muscle_D22_BP_set2a_0.05, n_top = 5, ncol = 2, order_levels = GeneSet_PairSet_target$set2)
GSEA_PublishDotplot(GSEA_GO_muscle_D22_BP_set2b_0.05, n_top = 5, ncol = 2, order_levels = GeneSet_PairSet_target$set2)
GSEA_PublishDotplot(GSEA_GO_muscle_D22_BP_set3a_0.05, n_top = 5, ncol = 2, order_levels = GeneSet_PairSet_target$set3)

GSEA_GO_muscle_D22_BP_set1_final_0.05 <- GSEA_GO_muscle_D22[GSEA_GO_muscle_D22$ONTOLOGY == "BP", ] %>%
  filter(contrast %in% c(
    "TM_N vs Con_N","Con_AHS vs Con_N",
    "TM_AHS vs Con_N","TM_AHS vs Con_AHS"
  ))
GSEA_PublishDotplot(GSEA_GO_muscle_D22_BP_set1_final_0.05, n_top = 5, ncol = 2, order_levels = GeneSet_PairSet_target$set1)

# ViSEAGO functional enrichment clustering - HS analysis -----------------------------------------------------------------

## define parameters

node_size <- 10

GeneSet_PairSet_target_v2 <- list(
  set1 = c(
    "TM_N vs Con_N",
    "Con_AHS vs Con_N",
    "Con_AHS vs TM_N",
    "TM_AHS vs Con_N",
    "TM_AHS vs TM_N",
    "TM_AHS vs Con_AHS"),
  set2 = c(
    "TM_N vs Con_N",
    "Con_CHS vs Con_N",
    "Con_CHS vs TM_N",
    "TM_CHS vs Con_N",
    "TM_CHS vs TM_N",
    "TM_CHS vs Con_CHS"))

## GO annotations from Bioconductor::org.Gg.eg.db

Bioc         <- ViSEAGO::Bioconductor2GO()
GENE2GO_map  <- ViSEAGO::annotate("org.Gg.eg.db", Bioc)

## Background unique non-NA genes (ENTREZ)

background_genes_muscle_D22_entrez_clean <- unique(na.omit(background_genes_muscle_D22_entrez))

## ViSEAGO cluster analysis for set 1 - Acute HS

pairwise_contrast_set1 <- GeneSet_PairSet_target_v2$set1

degs_sets_muscle_D22_entrez_set1 <- setNames(lapply(pairwise_contrast_set1, function(lbl){
  gs_ens <- DEG_full_info_muscle_D22 %>% filter(contrast == lbl) %>% pull(Gene)
  unique(na.omit(ens2entrez_map_tbl_mod[gs_ens]))
}), pairwise_contrast_set1)

ORA_topGO_muscle_D22_list_set1 <- list()

for (lbl in names(degs_sets_muscle_D22_entrez_set1)) {
  lbl_clean <- gsub("\\s+", "_", lbl)
  obj_go    <- paste0("BP_", lbl_clean)
  obj_res   <- paste0("elim_BP_", lbl_clean)
  
  gene_sel <- intersect(degs_sets_muscle_D22_entrez_set1[[lbl]], background_genes_muscle_D22_entrez_clean)
  if (length(gene_sel) == 0L) next
  
  assign(
    obj_go,
    ViSEAGO::create_topGOdata(
      geneSel  = gene_sel,
      allGenes = background_genes_muscle_D22_entrez_clean,
      gene2GO  = GENE2GO_map,
      ont      = "BP",
      nodeSize = node_size
    )
  )
  
  assign(
    obj_res,
    topGO::runTest(get(obj_go), algorithm = "elim", statistic = "fisher")
  )
  
  ORA_topGO_muscle_D22_list_set1[[lbl_clean]] <- c(obj_go, obj_res)
}

## GO terms semantic similarity

ViSEAGO_merge_muscle_D22_set1 <- ViSEAGO::merge_enrich_terms(Input = ORA_topGO_muscle_D22_list_set1)
ViSEAGO_goSS1_muscle_D22_set1 <- ViSEAGO::build_GO_SS(gene2GO = GENE2GO_map, enrich_GO_terms = ViSEAGO_merge_muscle_D22_set1)
ViSEAGO_goSS2_muscle_D22_set1 <- ViSEAGO::compute_SS_distances(ViSEAGO_goSS1_muscle_D22_set1, distance = "Wang")

ViSEAGO_hm_WangClu_WardD2_muscle_D22_set1 <- ViSEAGO::GOterms_heatmap(
  ViSEAGO_goSS2_muscle_D22_set1,
  showIC = TRUE,
  showGOlabels = TRUE,
  GO.tree = list(
    tree = list(distance = "Wang", aggreg.method = "ward.D2"),
    cut  = list(dynamic = list(pamStage = TRUE,
                               pamRespectsDendro = TRUE,
                               deepSplit = 2,
                               minClusterSize = 2))),
  samples.tree = NULL)

ViSEAGO::show_heatmap(ViSEAGO_hm_WangClu_WardD2_muscle_D22_set1, "GOterms")
ViSEAGO::show_table(ViSEAGO_merge_muscle_D22_set1)

## ViSEAGO cluster analysis for set 2 - Chronic HS

DEG_full_info_muscle_D22$contrast <- droplevels((DEG_full_info_muscle_D22$contrast))

pairwise_contrast_set2 <- GeneSet_PairSet_target_v2$set2

degs_sets_muscle_D22_entrez_set2 <- setNames(lapply(pairwise_contrast_set2, function(lbl){
  gs_ens <- DEG_full_info_muscle_D22 %>% filter(contrast == lbl) %>% pull(Gene)
  unique(na.omit(ens2entrez_map_tbl_mod[gs_ens]))
}), pairwise_contrast_set2)

ORA_topGO_muscle_D22_list_set2 <- list()

for (lbl in names(degs_sets_muscle_D22_entrez_set2)) {
  lbl_clean <- gsub("\\s+", "_", lbl)
  obj_go    <- paste0("BP_", lbl_clean)
  obj_res   <- paste0("elim_BP_", lbl_clean)
  
  gene_sel <- intersect(degs_sets_muscle_D22_entrez_set2[[lbl]], background_genes_muscle_D22_entrez_clean)
  if (length(gene_sel) == 0L) next
  
  assign(
    obj_go,
    ViSEAGO::create_topGOdata(
      geneSel  = gene_sel,
      allGenes = background_genes_muscle_D22_entrez_clean,
      gene2GO  = GENE2GO_map,
      ont      = "BP",
      nodeSize = node_size
    )
  )
  
  assign(
    obj_res,
    topGO::runTest(get(obj_go), algorithm = "elim", statistic = "fisher")
  )
  
  ORA_topGO_muscle_D22_list_set2[[lbl_clean]] <- c(obj_go, obj_res)
}

## GO terms semantic similarity

ViSEAGO_merge_muscle_D22_set2 <- ViSEAGO::merge_enrich_terms(Input = ORA_topGO_muscle_D22_list_set2)
ViSEAGO_goSS1_muscle_D22_set2 <- ViSEAGO::build_GO_SS(gene2GO = GENE2GO_map, enrich_GO_terms = ViSEAGO_merge_muscle_D22_set2)
ViSEAGO_goSS2_muscle_D22_set2 <- ViSEAGO::compute_SS_distances(ViSEAGO_goSS1_muscle_D22_set2, distance = "Wang")

ViSEAGO_hm_WangClu_WardD2_muscle_D22_set2 <- ViSEAGO::GOterms_heatmap(
  ViSEAGO_goSS2_muscle_D22_set2,
  showIC = TRUE,
  showGOlabels = TRUE,
  GO.tree = list(
    tree = list(distance = "Wang", aggreg.method = "ward.D2"),
    cut  = list(dynamic = list(pamStage = TRUE,
                               pamRespectsDendro = TRUE,
                               deepSplit = 2,
                               minClusterSize = 2))),
  samples.tree = NULL)

ViSEAGO::show_heatmap(ViSEAGO_hm_WangClu_WardD2_muscle_D22_set2, "GOterms")
ViSEAGO::show_table(ViSEAGO_merge_muscle_D22_set2)

# Supp_tables -------------------------------------------------------------

## convert all Supp_tables to proper data.frame class, if needed

TxiSalmon_muscle_abund <- as.data.frame(TxiSalmon_muscle$abundance) %>%
  rownames_to_column("isoform_id") %>%
  rename_with(~ paste0(.x, "_TPM"), -isoform_id)
TxiSalmon_muscle_counts <- as.data.frame(TxiSalmon_muscle$counts) %>%
  rownames_to_column("isoform_id") %>%
  rename_with(~ paste0(.x, "_count"), -isoform_id)
TxiSalmon_muscle_ExpMatrix <- inner_join(TxiSalmon_muscle_abund, TxiSalmon_muscle_counts, by = "isoform_id")
rm(TxiSalmon_muscle_abund,TxiSalmon_muscle_counts)

ViSEAGO_res_set1_df <- as.data.frame(ViSEAGO_merge_muscle_D22_set1@data)
ViSEAGO_res_set2_df <- as.data.frame(ViSEAGO_merge_muscle_D22_set2@data)

TxiSalmon_muscle_ExpMatrix                <- as.data.frame(TxiSalmon_muscle_ExpMatrix)
DEG_full_info_muscle_D22_all_anno         <- as.data.frame(DEG_full_info_muscle_D22_all_anno)
ORA_GO_muscle_D22                         <- as.data.frame(ORA_GO_muscle_D22)
GSEA_GO_muscle_D22                        <- as.data.frame(GSEA_GO_muscle_D22)
ViSEAGO_res_set1_df                       <- as.data.frame(ViSEAGO_res_set1_df)
ViSEAGO_res_set2_df                       <- as.data.frame(ViSEAGO_res_set2_df)

## create the excel file table of content

Rdataframes_titles <- data.frame(
  Sheet = c(
    "Sheet1",
    "Sheet2",
    "Sheet3",
    "Sheet4",
    "Sheet5",
    "Sheet6"
  ),
  
  DataFrameName = c(
    "TxiSalmon_muscle_ExpMatrix",
    "DEG_full_info_muscle_D22_all_anno",
    "ORA_GO_muscle_D22",
    "GSEA_GO_muscle_D22",
    "ViSEAGO_res_set1_df",
    "ViSEAGO_res_set2_df"
  ),
  
  DataFrameNote = c(
    "Salmon expression matrix results at the level of transcripts with raw and TPM counts",
    "Differentially Expressed Genes (DEGs) for 6 group comparisons - Day 22 post-hatch",
    "ORA for GO terms results (Biological Process) - 6 group comparisons - Day 22 post-hatch",
    "GSEA for GO terms results (Biological Process) - 6 group comparisons - Day 22 post-hatch",
    "ViSEAGO GO BP clsuter results for acute heat stress groups",
    "ViSEAGO GO BP clsuter results for chronic heat stress groups"
  )
)

## create the excel file

write_xlsx(list(
  
  ContentTable = Rdataframes_titles,
  
  Sheet1  = TxiSalmon_muscle_ExpMatrix,
  Sheet2  = DEG_full_info_muscle_D22_all_anno,
  Sheet3  = ORA_GO_muscle_D22,
  Sheet4  = GSEA_GO_muscle_D22,
  Sheet5  = ViSEAGO_res_set1_df,
  Sheet6  = ViSEAGO_res_set2_df
),

file.path(sys_dir,"TMBroilers_Transcriptomics/muscle/supp_tables/SuppTables_TMBroilersTranscriptomics_Muscle.xlsx")
)


