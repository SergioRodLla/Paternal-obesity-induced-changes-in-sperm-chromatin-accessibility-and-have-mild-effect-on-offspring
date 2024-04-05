######################################################################################################################################
## This R code was used to perform peak differential accessibility analysis and prepare inputs for HOMER TF motif analysis for
## INSERT CITATION + DOI
## DISCLAIMER: Manuscript figures may differ from the ones obtained with this code ONLY for enhancing visualization of final figures.
## Gene Ontology enrichment bubble plots were filtered to show non-redundant informative BP pathways.
## Sergio R. Llana 5/04/2024
######################################################################################################################################

# IMPORTING PACKAGES ####################
library(DESeq2)
library(ggplot2)
library(plotly)
library(tidyverse)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(gprofiler2)
library(clusterProfiler)
library(forcats)
library(circlize)

# IMPORTING PEAK TABLES ####################

# 0/1 matrix indicating presence of peaks
overs <- read.table("consensus_overlaps.tsv",
                    header = T)
# Peak count matrix
counts <- read.table("consensus_counts.tsv",
                     header = T)
# Keep peaks that appear in at least two samples
counts_clean <- counts[which(rowSums(overs[,4:ncol(overs)])>=2),]
rownames(counts_clean) <- paste0(counts_clean$chr,":",counts_clean$start,"-",counts_clean$end)
colnames(counts_clean) <- gsub(".THSS.bam", "", colnames(counts_clean))

# Load metadata file: Table containing information about each sample (e.g. diet)
enh <- readxl::read_xlsx("metadata.xlsx") %>%
  dplyr::arrange(ID) %>%
  select(-ID)

rownames(enh) <- colnames(counts_clean[,4:ncol(counts_clean)])

# RUN PCA ####################

dds <- DESeqDataSetFromMatrix(counts_clean[,4:ncol(counts_clean)], enh, ~diet)

countsraw <- counts(dds)
countsvst <- vst(countsraw)

# PCA
pc <- prcomp(t(countsvst))
pcx <- as.data.frame(pc$x)
#head(summary(pc)$importance[2,])

ggplot(pcx,aes(pcx[,1],pcx[,2],colour=as.factor(enh$diet),label=rownames(enh))) +
  geom_point(size=4) +
  labs(x=paste0("PC1 ", round((summary(pc)$importance[2,1])*100,2)),y=paste0("PC2 ", round((summary(pc)$importance[2,2])*100,2)),
       colour="diet") + 
  theme_minimal() + 
  geom_text(vjust=1.5)

# Observe PC individually
pca.res.df <- pc$x[,1:5] %>%
  as_tibble() %>%
  add_column(sample = rownames(enh),
             diet = enh$diet)

pca.pivot <- pivot_longer(pca.res.df,
                          cols = PC1:PC5,
                          names_to = "PC",
                          values_to = "loadings")

pca.pivot$PC <- factor(pca.pivot$PC, levels = c("PC1", "PC2", "PC3", "PC4", "PC5"))

PC_num <- c("PC1", "PC2", "PC3", "PC4", "PC5")
PC_var <- as.character(round((summary(pc)$importance[2,1:5])*100,2))
PC_tags <- paste0(PC_num, " : ", PC_var)
names(PC_tags) <- c("PC1", "PC2", "PC3", "PC4", "PC5")

ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill=diet) + 
  geom_bar(stat="identity") +
  facet_wrap(~PC, labeller = labeller(PC = PC_tags)) +
  labs(title="Individual PCs") +
  theme_bw() +
  coord_flip()

# RUN PEAK DIFFERENTIAL ACCESSIBILITY ANALYSIS ####################

keep <- rowSums(counts(dds)) >= 4 #Set threshold for peak counts
dds_clean <- dds[keep,]
dds_clean <- DESeq(dds_clean) # Run DESeq2

#resultsNames(dds_clean)

# Diff expr WD vs CD
res <- results(dds_clean, contrast=c("diet", "WD", "CD"))

# Annotate peaks to genomic regions using TSS
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
counts_clean_gr <- makeGRangesFromDataFrame(counts_clean)
peakAnno_all <- annotatePeak(counts_clean_gr, tssRegion=c(-3000, 3000),
                             TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno_all.df <- as.data.frame(peakAnno_all) %>%
  mutate(peakID = paste0(seqnames, ":", start, "-", end))

res.df <- as_tibble(res, rownames = "peaks") %>%
  dplyr::mutate(signif = ifelse(pvalue<0.05& log2FoldChange < -1, "Closed", ifelse(pvalue<0.05&log2FoldChange > 1, "Open", "NS" )))

res.df <- left_join(res.df, peakAnno_all.df, by = join_by(peaks == peakID))

# Run to get peaks with significant increased/decreased accessibility 
#View(res.df %>% filter(signif %in% c("Open", "Closed")))

# VOLCANO PLOT ###################

signif.colors <- c(Closed = "#2C467A", Open = "#BE684D", NS = "#000000")

#theme_old <- theme_set(theme_classic())

vplot <- ggplot(res.df) +
  aes(y=-log10(pvalue), x=log2FoldChange, text = paste("Symbol:", SYMBOL), color = signif) +
  geom_point(size=2) +
  scale_color_manual(values = signif.colors, guide="none") +
  scale_x_continuous(limits = c(-4.5, 4.5)) +
  guides(color=guide_legend(title="Diff. Accessibility")) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", linewidth=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  annotate("text", x = 3, y = 0, label=paste0("Open = ", nrow(res.df %>% filter(signif=="Open"))), color="#BE684D") +
  annotate("text", x = -3, y = 0, label=paste0("Closed = ", nrow(res.df %>% filter(signif=="Closed"))), color="#2C467A") +
  xlab(bquote(log[2]*"FC")) +
  ylab(bquote(-log[10]*"P.value"))

vplot

# Run interactively
ggplotly(vplot)

# GENERATE PEAK BED FILES FOR TF MOTIF ENRICHMENT ###################

# Get peaks with increased accessibility (both WD and CD conditions) for HOMER tool

# Get open/closed peaks, using thresholds
open_peaks.df <- res.df %>%
  filter(!grepl("chrG|chrJ", seqnames)) %>%
  filter(signif == "Open") %>%
  mutate(annotation_clean = gsub("\\s*\\([^\\)]+\\)", "", annotation)) %>%
  mutate(annotation_clean = gsub(" ", "", annotation_clean)) %>%
  mutate(annotation_clean = gsub("'", "", annotation_clean)) %>%
  # Create a unique peak ID (for HOMER tool)
  mutate(Unique_PID = paste0(SYMBOL,"_",annotation_clean)) %>%
  select(seqnames, start, end, Unique_PID) %>%
  mutate(fillScore=0) %>%
  mutate(fillStrand=".") %>%
  na.omit()

closed_peaks.df <- res.df %>%
  filter(!grepl("chrG|chrJ", seqnames)) %>%
  filter(signif == "Closed") %>%
  mutate(annotation_clean = gsub("\\s*\\([^\\)]+\\)", "", annotation)) %>%
  mutate(annotation_clean = gsub(" ", "", annotation_clean)) %>%
  mutate(annotation_clean = gsub("'", "", annotation_clean)) %>%
  mutate(Unique_PID = paste0(SYMBOL,"_",annotation_clean)) %>%
  select(seqnames, start, end, Unique_PID) %>%
  mutate(fillScore=0) %>%
  mutate(fillStrand=".") %>%
  na.omit()

# Export using BED-like format
write.table(open_peaks.df, "signif_open_peaks.bed", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(closed_peaks.df, "signif_closed_peaks.bed", quote = F, col.names = F, row.names = F, sep = "\t")

# PLOT PEAKS PER GENOMIC REGION ###################

res_clean_signif.df <- res.df %>%
  mutate(annotation_clean = gsub("\\s*\\([^\\)]+\\)", "", annotation)) %>%
  filter(signif %in% c("Open", "Closed"))

res_clean_signif_pie.df <- res_clean_signif.df %>%
  arrange(desc(annotation_clean)) %>%
  group_by(signif, annotation_clean) %>%
  summarize(count = n()) %>%
  group_by(signif) %>%
  mutate(n_group = sum(count)) %>%
  mutate(perc = round(count/n_group*100,2))

res_clean_signif.df <- left_join(res_clean_signif.df, res_clean_signif_pie.df, by = "annotation_clean")

#theme_old <- theme_set(theme_classic())
signif.colors <- c(Closed = "#2C467A", Open = "#BE684D", NS = "#000000")

# Barplot
ggplot(res_clean_signif_pie.df, aes(x=annotation_clean, y=count, fill=signif)) + 
  geom_bar(stat="identity", position = position_dodge(), color="black") + 
  coord_flip() + scale_fill_discrete(name="Diff. Accessibility") +
  scale_fill_manual(values = signif.colors)

# Pie chart
ggplot(res_clean_signif_pie.df, aes(x="", y=perc, fill=annotation_clean)) +
  geom_bar(stat="identity", width=1, color="white") +
  facet_wrap(~signif) +
  coord_polar("y", start=0) +
  theme_void() +
  geom_text(aes(label = paste0(round(perc, 2),"%")), position = position_stack(vjust = 0.5), size=3) +
  scale_fill_brewer(palette="Set1", name = "Genomic regions")


# GENE ONTOLOGY ENRICHMENT ###################

# get the gene names of diff. open and closed peaks
diff_genes_open <- res.df %>%
  filter(signif == "Open") %>%
  select(SYMBOL)
diff_genes_closed <- res.df %>%
  filter(signif == "Closed") %>%
  select(SYMBOL)

# Run g:Profiler g:GOSt for genes associated with OPEN REGIONS
gprofiler_open <- gost(unique(diff_genes_open$SYMBOL), organism = "mmusculus", correction_method = "fdr", evcodes = T)
gostplot(gprofiler_open)

# Filter for Biological Process
gprofiler_open_bp <- gprofiler_open$result %>%
  filter(source == "GO:BP")

# Subset pathways to remove GO terms redundancy (related to GO hierarchy)
nervous_system <- c("sensory neuron axon guidance", "synapse organization", "neuron differentiation", "neurogenesis", "cell projection organization")
digestive_tract <- c("embryonic digestive tract development", "embryonic digestive tract morphogenesis", "cell junction organization", "cell differentiation")

# Create bubble plot
gprofiler_open_bp %>%
  filter(term_name %in% c(digestive_tract, nervous_system)) %>%
  mutate(term_name = factor(term_name, rev(c(digestive_tract, nervous_system)))) %>%
  mutate(Category = ifelse(term_name %in% digestive_tract, "Digestive system", "Nervous system")) %>%
  mutate(Category = factor(Category, c("Digestive system", "Nervous system"))) %>%
  ggplot(., aes(x=p_value, y=term_name, size=intersection_size/term_size)) +
  geom_point(aes(color=Category)) +
  scale_color_manual(values=c("#FF61CC", "#7CAE00")) +
  guides(size=guide_legend(title="Gene coverage")) +
  labs(x="adj.P.value (FDR)", y="Biological Process") +
  scale_x_reverse() +
  coord_cartesian(xlim = c(0.03, 0))

# Run g:Profiler g:GOSt for genes associated with CLOSED REGIONS
gprofiler_closed <- gost(unique(diff_genes_closed$SYMBOL), organism = "mmusculus", correction_method = "fdr", evcodes = T)
gostplot(gprofiler_closed)

gprofiler_closed_bp <- gprofiler_closed$result %>%
  filter(source == "GO:BP")

heart_cardio <- c("regulation of secondary heart field cardioblast proliferation", "heart contraction", "cardiac cell fate commitment", "regulation of angiogenesis", "cardioblast differentiation", "cell communication involved in cardiac conduction")
transcription <- c("peptidyl-serine phosphorylation", "canonical NF-kappaB signal transduction", "RNA biosynthetic process", "DNA-templated transcription", "transcription by RNA polymerase II", "protein phosphorylation")
nervous_system <- c("regulation of timing of neuron differentiation", "forebrain generation of neurons", "neuron projection development", "forebrain radial glial cell differentiation", "calcium ion transport")
extras <- c("tube development", "regulation of fat cell differentiation", "glucose transmembrane transport", "regulation of insulin secretion", "regulation of primary metabolic process")

gprofiler_closed_bp %>%
  filter(term_name %in% c(heart_cardio, nervous_system, transcription, extras)) %>%
  mutate(term_name = factor(term_name, rev(c(heart_cardio, nervous_system, transcription, extras)))) %>%
  mutate(Category = ifelse(term_name %in% heart_cardio, "Cardiovascular devel.",
                           ifelse(term_name %in% nervous_system, "Nervous system",
                                  ifelse(term_name %in% transcription, "Transcription", "Other")))) %>%
  mutate(Category = factor(Category, c("Cardiovascular devel.", "Nervous system", "Transcription", "Other"))) %>%
  ggplot(., aes(x=p_value, y=term_name, size=intersection_size/term_size)) +
  geom_point(aes(color=Category)) +
  guides(size=guide_legend(title="Gene coverage")) +
  labs(x="adj.P.value (FDR)", y="Biological Process") +
  scale_x_reverse() +
  coord_cartesian(xlim = c(0.05, 0))

# TF MOTIFS CIRCLE PLOTS ###################

### Cirlcle plot for CLOSED regions' TF motifs

# Pathways were inferred from GO annotations from Mouse Genome Informatics database. TFs were grouped accordingly.
# Notice that Nfe2l2 corresponds to TF morif "ARE(NR)" and PPAR to TF motif "PPARE(NR),DR1"
Cardio_dev <- c("ERG", "Foxh1", "Oct4:Sox17", "TR4", "Nfe2l2")
Nerv_dev <- c("CUX1", "ETV1", "ETV4", "TR4", "bHLHE40", "Nfe2l2")
Skeletal_muscle_dev <- c("Nur77", "Sox15")
Fat_tissue_dev <- c("Nur77", "PPAR")
Metab <- c("PPAR", "Nfe2l2")

TF_closed.list <- unique(c(Cardio_dev, Nerv_dev, Skeletal_muscle_dev, Fat_tissue_dev, Metab))

group_list <- list("Cardiovascular development" = Cardio_dev,
                   "Nervous system development" = Nerv_dev, 
                   "Skeletal muscle development" = Skeletal_muscle_dev,
                   "Fat tissue development" = Fat_tissue_dev,
                   "Metabolism" = Metab)

TF_closed_grouped <- matrix(0, nrow = length(group_list), ncol = length(TF_closed.list), dimnames = list(names(group_list), names(TF_closed.list)))

for (i in seq_along(group_list)) {
  TF_closed_grouped[i, group_list[[i]]] <- 1
}
TF_closed_grouped.df <- TF_closed_grouped %>% as.data.frame()


#set q-value for the TFs (from HOMER results)
qvals <- c("bHLHE40"=0.0303, "CUX1"=0.0303, "ERG"=0.0303, "ETV1"=0.0211, "ETV4"=0.0303, "Foxh1"=0.0303, "Nur77"=0.0303, "Nfe2l2"=0.0303,"PPAR"=0.0303,"Sox15"=0.0437, "Oct4:Sox17"=0.0016, "TR4"=0.0472)

#define color range 
## Use n equally spaced breaks to assign each value to n-1 equal sized bins 
ii <- cut(qvals, breaks = seq(min(qvals), 0.05, len = 100), 
          include.lowest = TRUE)
## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
colors <- colorRampPalette(c("darkblue", "#b9f8ff"))(99)[ii]
names(colors) <- c("bHLHE40", "CUX1" , "ERG" , "ETV1", "ETV4", "Foxh1", "Nur77", "Nfe2l2", "PPAR", "Sox15", "Oct4:Sox17", "TR4")

# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), width=1.75, height=5, title='') {
  scale = (length(lut)-1)/(max-min)
  
  #dev.new(width, height)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

# Assign colors to grids
grid.col <- c("Cardiovascular development" = "#F8766D", 
              "Nervous system development" = "#A3A500", 
              "Skeletal muscle development" = "#00BF7D",
              "Fat tissue development" = "#00B0F6",
              "Metabolism" = "#E76BF3",
              "bHLHE40" = unname(colors["bHLHE40"]), 
              "CUX1" = unname(colors["CUX1"]), 
              "ERG" = unname(colors["ERG"]), 
              "ETV1" = unname(colors["ETV1"]), 
              "ETV4" = unname(colors["ETV4"]), 
              "Foxh1" = unname(colors["Foxh1"]), 
              "Nur77" = unname(colors["Nur77"]), 
              "Nfe2l2" = unname(colors["Nfe2l2"]), 
              "PPAR" = unname(colors["PPAR"]),  
              "Sox15"= unname(colors["Sox15"]), 
              "Oct4:Sox17" = unname(colors["Oct4:Sox17"]), 
              "TR4" = unname(colors["TR4"]))

# Generate circle plot
chordDiagramFromMatrix(TF_closed_grouped, annotationTrack = c("grid"), grid.col = grid.col, transparency = 0.3, directional = -1, preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(TF_closed_grouped))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

# Generate color bar to indicate motif q-value
color.bar(colorRampPalette(c("darkblue", "#b9f8ff"))(99), min=0, max=0.05, nticks = 2, title="q-value", width = 0.5, height = 3)



### Cirlcle plot for OPEN regions' TF motifs

# Note that Jun corresponds to TF motif "Jun-AP1(bZIP)"
Dev_pro <- c("Egr2", "Jun", "NFIL3", "Npas4", "RUNX1", "Srebp1a", "Srebp2", "Tcfcp2l1", "Zic2", "Zic3", "ZNF143")
Sig_stim <- c("Egr2", "IRF3", "Jun", "KLF14", "Npas4", "RUNX1", "Srebp1a")
Metab <- c("Egr2", "Srebp1a", "Srebp2")

TF_open.list <- unique(c(Dev_pro, Sig_stim, Metab))

group_names <- c("Developmental process", "Signaling", "Metabolism")

group_list <- list("Developmental process" = Dev_pro, "Signaling" = Sig_stim, "Metabolism" = Metab)

TF_open_grouped <- matrix(0, nrow = length(group_list), ncol = length(TF_open.list), dimnames = list(names(group_list), names(TF_open.list)))

for (i in seq_along(group_list)) {
  TF_open_grouped[i, group_list[[i]]] <- 1
}
TF_open_grouped.df <- TF_open_grouped %>% as.data.frame()


#set q-value for the TFs
qvals <- c("Egr2" = 0.0084, "IRF3" = 0.0176, "Jun" = 0.0190, "KLF14" = 0.0190, "NFIL3" = 0.0190, "Npas4" = 0.0089, "RUNX1" = 0.0112, "Srebp1a" = 0.0224, "Srebp2" = 	0.0190,  "Tcfcp2l1"= 	0.0084, "Zic2" = 0.0081, "Zic3" = 0.0000, "ZNF143" = 0.0042)

ii <- cut(qvals, breaks = seq(min(qvals), 0.03, len = 100), # set max = 0.03 to better distinguish q-val differences
          include.lowest = TRUE)
colors <- colorRampPalette(c("darkblue", "#b9f8ff"))(99)[ii]
names(colors) <- c("Egr2", "IRF3" , "Jun" , "KLF14", "NFIL3", "Npas4", "RUNX1", "Srebp1a", "Srebp2", "Tcfcp2l1", "Zic2", "Zic3", "ZNF143")

grid.col <- c("Developmental process" = "#F8766D", "Signaling" = "#7CAE00", "Metabolism" = "#00BFC4", 
              "Egr2" = unname(colors["Egr2"]), 
              "IRF3" = unname(colors["IRF3"]), 
              "Jun" = unname(colors["Jun"]), 
              "KLF14" = unname(colors["KLF14"]), 
              "NFIL3" = unname(colors["NFIL3"]), 
              "Npas4" = unname(colors["Npas4"]), 
              "RUNX1" = unname(colors["RUNX1"]), 
              "Srebp1a" = unname(colors["Srebp1a"]), 
              "Srebp2" = unname(colors["Srebp2"]),  
              "Tcfcp2l1"= unname(colors["Tcfcp2l1"]), 
              "Zic2" = unname(colors["Zic2"]), 
              "Zic3" = unname(colors["Zic3"]), 
              "ZNF143" = unname(colors["ZNF143"]))


chordDiagramFromMatrix(TF_open_grouped, annotationTrack = c("grid"), grid.col = grid.col, transparency = 0.3, directional = -1, preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(TF_open_grouped))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

color.bar(colorRampPalette(c("darkblue", "#b9f8ff"))(99), min=0, max=0.03, nticks = 0, title="q-value", width = 0.5, height = 3)






