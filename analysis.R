# carga de librerías
require(DESeq2)
require(sva)
require(ggplot2)
require(VennDiagram)
require(org.Hs.eg.db)
require(biomaRt)
require(clusterProfiler)
require(reshape)
require(RColorBrewer)
require(pheatmap)

# 1. LECTURA Y PRE-PROCESADO DE LOS DATOS DE RNA-SEQ
covariables <- read.csv2('./data/targets.csv', sep = ',')
covariables$Group <- factor(covariables$Group, levels = c('NIT', 'SFI', 'ELI'))
covariables$sex <- factor(covariables$sex, levels = c('male', 'female'))
covariables$Sample_Name <- gsub('-', '.', covariables$Sample_Name)
contajes <- read.csv('./data/counts.csv', sep = ';', header = TRUE, row.names = 1)
# se escogen 10 muestras aleatorias de cada tipo de infiltración en el tiroides
set.seed(1234)
muestras_NIT <- sample(which(covariables$Group == 'NIT'), size = 10, replace = FALSE)
muestras_SFI <- sample(which(covariables$Group == 'SFI'), size = 10, replace = FALSE)
muestras_ELI <- sample(which(covariables$Group == 'ELI'), size = 10, replace = FALSE)
covariables_muestra <- covariables[c(muestras_NIT, muestras_SFI, muestras_ELI), ]
contajes_muestra <-  contajes[, covariables_muestra$Sample_Name]
rownames(contajes_muestra) <- gsub('\\..*', '', rownames(contajes_muestra))


# 2. INSPECCIÓN DE LOS DATOS
# dado que las distribuciones de los contajes son fuertemente asimétricas (debido a pocos contajes muy elevados), 
# se utiliza la transformación log2 para aproximadamente normalizar las mismas
pseudo_contajes_muestra <- log2(contajes_muestra + 1)

# boxplots
df <- melt(pseudo_contajes_muestra, variable_name = 'muestra')
df <- merge(df, covariables_muestra[, c('Sample_Name', 'Group')], by.x = 'muestra', by.y = 'Sample_Name', all = TRUE)
colnames(df)[3] <- 'grupo'
ggplot(df, aes(x = muestra, y = value, fill = grupo)) + geom_boxplot() + 
  xlab('') + ylab(expression(log[2](recuento + 1))) + 
  scale_fill_manual(values = c('#274948', '#a57b82', '#e39e2f')) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# gráficos de densidad
ggplot(df, aes(x = value, colour = muestra, fill = muestra)) + 
  geom_density(alpha = 0.2, size = 1.25) + facet_wrap(~ grupo) +
  theme(legend.position = 'none') + xlab(expression(log[2](recuento + 1)))

# PCA
dds <- DESeqDataSetFromMatrix(countData = contajes_muestra, colData = covariables_muestra, design = ~ Group + sex)
dds <- dds[rowSums(counts(dds)) > 1, ] # filtraje
vsd <- vst(dds, blind = FALSE) # para hacer los recuentos aproximadamente homocedásticos
plotPCA(vsd, intgroup = c('Group', 'sex')) + 
  scale_color_manual(values = c('#e39e2f', '#fdde9b', '#274948', '#cce5e5', '#a57b82', '#ffe9eb')) +
  theme_classic()

# heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Group)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)


# 3. NORMALIZACIÓN
# Relative Log Expression (RLE)
dds <- estimateSizeFactors(dds)


# 4. ANÁLISIS DE EXPRESIÓN DIFERENCIAL
# estimación de las dispersiones
dds <- estimateDispersions(dds)
# cálculo de la expresión diferencial haciendo uso de la distribución discreta binomial negativa
dds <- nbinomWaldTest(dds)
# control de efectos batch
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ Group + sex, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv = 1)
dds$SV1 <- svseq$sv[,1]
design(dds) <- ~ SV1 + Group + sex

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

res_ELI_vs_NIT <- results(dds, contrast = c('Group', 'ELI', 'NIT'))
res_ELI_vs_NIT_sig <- subset(res_ELI_vs_NIT, padj < 0.01 & abs(log2FoldChange) > 1.5)
res_ELI_vs_SFI <- results(dds, contrast = c('Group', 'ELI', 'SFI'))
res_ELI_vs_SFI_sig <- subset(res_ELI_vs_SFI, padj < 0.01 & abs(log2FoldChange) > 1.5)
res_SFI_vs_NIT <- results(dds, contrast = c('Group', 'SFI', 'NIT'))
res_SFI_vs_NIT_sig <- subset(res_SFI_vs_NIT, padj < 0.01 & abs(log2FoldChange) > 1.5)


# ANOTACIÓN DE RESULTADOS
mart <- useDataset('hsapiens_gene_ensembl', useMart('ENSEMBL_MART_ENSEMBL'))

getSymbols <- function(ensemblIDs) {
  annotLookup <- getBM(
    mart = mart,
    attributes = c('ensembl_gene_id', 'external_gene_name'),
    filter = 'ensembl_gene_id',
    values = ensemblIDs,
    uniqueRows = TRUE)
  return(annotLookup)
}

res_ELI_vs_NIT_sig$entrez <- mapIds(org.Hs.eg.db,
                                    keys=row.names(res_ELI_vs_NIT_sig),
                                    column='ENTREZID',
                                    keytype='ENSEMBL',
                                    multiVals='first')
res_ELI_vs_NIT_sig <- merge(data.frame(res_ELI_vs_NIT_sig), getSymbols(rownames(res_ELI_vs_NIT_sig)), by.x = 0, by.y = 'ensembl_gene_id', all = TRUE)
colnames(res_ELI_vs_NIT_sig)[c(1,9)] <- c('ensembl', 'symbol')
head(res_ELI_vs_NIT_sig[order(abs(res_ELI_vs_NIT_sig$log2FoldChange), decreasing = TRUE),], 20)

res_ELI_vs_SFI_sig$entrez <- mapIds(org.Hs.eg.db,
                                    keys=row.names(res_ELI_vs_SFI_sig),
                                    column='ENTREZID',
                                    keytype='ENSEMBL',
                                    multiVals='first')
res_ELI_vs_SFI_sig <- merge(data.frame(res_ELI_vs_SFI_sig), getSymbols(rownames(res_ELI_vs_SFI_sig)), by.x = 0, by.y = 'ensembl_gene_id', all = TRUE)
colnames(res_ELI_vs_SFI_sig)[c(1,9)] <- c('ensembl', 'symbol')
head(res_ELI_vs_SFI_sig[order(abs(res_ELI_vs_SFI_sig$log2FoldChange), decreasing = TRUE),], 20)

res_SFI_vs_NIT_sig$entrez <- mapIds(org.Hs.eg.db,
                                    keys=row.names(res_SFI_vs_NIT_sig),
                                    column='ENTREZID',
                                    keytype='ENSEMBL',
                                    multiVals='first')
res_SFI_vs_NIT_sig <- merge(data.frame(res_SFI_vs_NIT_sig), getSymbols(rownames(res_SFI_vs_NIT_sig)), by.x = 0, by.y = 'ensembl_gene_id', all = TRUE)
colnames(res_SFI_vs_NIT_sig)[c(1,9)] <- c('ensembl', 'symbol')
head(res_SFI_vs_NIT_sig[order(abs(res_SFI_vs_NIT_sig$log2FoldChange), decreasing = TRUE),], 20)


# 5. VISUALIZACIÓN DE RESULTADOS
# volcano plots
volcanoPlot <- function(resultados, lfc = 1.5, pval = 0.01) {
  tab <- data.frame(logFC = resultados$log2FoldChange, negLogPval = -log10(resultados$padj))
  plot(tab, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~adjustedpvalue), axes = FALSE)
  signGenes_up <- (tab$logFC > lfc & tab$negLogPval > -log10(pval))
  signGenes_down <- (tab$logFC < -lfc & tab$negLogPval > -log10(pval))
  points(tab[signGenes_up, ], pch = 16, cex = 0.8, col = '#d6aaaa')
  points(tab[signGenes_down, ], pch = 16, cex = 0.8, col = '#cce5e5')
  abline(h = -log10(pval), col = "gray", lty = 2)
  abline(v = c(-lfc, lfc), col = "gray", lty = 2)
  mtext(paste("FDR =", pval), side = 2, at = -log10(pval), cex = 0.6, line = 0.5, las = 1)
  mtext(c(expression(-1.5~log[2]~FC), expression(1.5~log[2]~FC)), side = 3, at = c(-3.5, 3.5), cex = 0.6, line = 0.5)
  axis(1, seq(-30, 30, 5), cex.axis = 0.7)
  axis(2, seq(0, 50, 10), cex.axis = 0.7)
}

volcanoPlot(res_ELI_vs_NIT)
volcanoPlot(res_ELI_vs_SFI)
volcanoPlot(res_SFI_vs_NIT)

# diagramas de Venn
venn.diagram(
  x = list(unique(res_ELI_vs_NIT_sig[res_ELI_vs_NIT_sig$log2FoldChange > 1.5,]$ensembl), 
           unique(res_ELI_vs_SFI_sig[res_ELI_vs_SFI_sig$log2FoldChange > 1.5,]$ensembl)),
  category.names = c('ELI vs NIT', 'ELI vs SFI'),
  filename = 'vennDiagram_ELIvsNIT_ELIvsSFI_up.png',
  output=TRUE,
  height = 430, 
  width = 430,
  resolution = 200,
  compression = 'lzw',
  lwd = 2,
  lty = 'blank',
  fill = c('#d6aaaa', '#ffdfe1'),
  cex = .6,
  fontface = 'bold',
  fontfamily = 'sans',
  cat.cex = 0.4,
  cat.fontface = 'bold',
  cat.default.pos = 'outer',
  cat.fontfamily = 'sans')

venn.diagram(
  x = list(unique(res_ELI_vs_NIT_sig[res_ELI_vs_NIT_sig$log2FoldChange < 1.5,]$ensembl), 
           unique(res_ELI_vs_SFI_sig[res_ELI_vs_SFI_sig$log2FoldChange < 1.5,]$ensembl)),
  category.names = c('ELI vs NIT', 'ELI vs SFI'),
  filename = 'vennDiagram_ELIvsNIT_ELIvsSFI_down.png',
  output=TRUE,
  height = 430, 
  width = 430,
  resolution = 150,
  compression = 'lzw',
  lwd = 2,
  lty = 'blank',
  fill = c('#274948', '#cce5e5'),
  cex = .8,
  fontface = 'bold',
  fontfamily = 'sans',
  cat.cex = 0.4,
  cat.fontface = 'bold',
  cat.default.pos = 'outer',
  cat.fontfamily = 'sans')

# análisis de significación biológica
# molecular function
res_ELI_vs_ALL_sig_upregulated_ego_MF <- enrichGO(
  gene = intersect(unique(res_ELI_vs_NIT_sig[res_ELI_vs_NIT_sig$log2FoldChange > 1.5,]$ensembl),
                   unique(res_ELI_vs_SFI_sig[res_ELI_vs_SFI_sig$log2FoldChange > 1.5,]$ensembl)),
  OrgDb = org.Hs.eg.db,
  ont = 'MF',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.01,
  readable = TRUE,
  keyType = 'ENSEMBL')

ggplot(data.frame(GO = res_ELI_vs_ALL_sig_upregulated_ego_MF$Description,
                  p.adjust = -log10(res_ELI_vs_ALL_sig_upregulated_ego_MF$p.adjust)),
       aes(x = p.adjust, y = reorder(GO, p.adjust))) +
  geom_bar(stat = 'identity', fill = alpha('#4c4c4c', 0.7)) + 
  xlab(expression(paste(-log[10], '(', italic('p-valor'), ' ajustado)'))) + ylab(NULL) +
  theme(plot.title = element_text(hjust = 0.5))
# biological process
res_ELI_vs_ALL_sig_upregulated_ego_BP <- enrichGO(
  gene = intersect(unique(res_ELI_vs_NIT_sig[res_ELI_vs_NIT_sig$log2FoldChange > 1.5,]$ensembl),
                   unique(res_ELI_vs_SFI_sig[res_ELI_vs_SFI_sig$log2FoldChange > 1.5,]$ensembl)),
  OrgDb = org.Hs.eg.db,
  ont = 'BP',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.01,
  readable = TRUE,
  keyType = 'ENSEMBL')

ggplot(data.frame(GO = res_ELI_vs_ALL_sig_upregulated_ego_BP$Description[1:50],
                  p.adjust = -log10(res_ELI_vs_ALL_sig_upregulated_ego_BP$p.adjust[1:50])),
       aes(x = p.adjust, y = reorder(GO, p.adjust))) +
  geom_bar(stat = 'identity', fill = alpha('#4c4c4c', 0.7)) + 
  xlab(expression(paste(-log[10], '(', italic('p-valor'), ' ajustado)'))) + ylab(NULL) +
  theme(plot.title = element_text(hjust = 0.5))
