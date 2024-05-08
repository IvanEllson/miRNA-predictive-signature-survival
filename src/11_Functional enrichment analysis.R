### Functional enrichment analysis

#### Load results from TAM2 (http://www.lirmed.com/tam2/) ####

Enrichment.TAM2_result <- read.delim("./results/Functional analysis TAM2/Results table.txt")

#----

#### Functional enrichment analysis of new signature ####

### Functions

functions.TAM2 <- Enrichment.TAM2_result[Enrichment.TAM2_result$Category=="Function",]
functions.TAM2 <- functions.TAM2[order(functions.TAM2$Term, decreasing = F),]
functions.TAM2 <- functions.TAM2[order(functions.TAM2$P.value, decreasing = F),]
functions.TAM2_top <- functions.TAM2[1:10,]


## Barplot
functions.TAM2_top$TermCount <- paste0(functions.TAM2_top$Term, " (", functions.TAM2_top$Count, ")")
functions.TAM2_top$logpval <- -log10(functions.TAM2_top$P.value)
functions.TAM2_top$TermCount <- factor(functions.TAM2_top$TermCount, levels = unique(functions.TAM2_top$TermCount))

plotbarF <- ggplot(functions.TAM2_top, aes(x=logpval, y=TermCount)) +
  geom_bar(stat = "identity") +
  scale_y_discrete(limits = rev(levels(functions.TAM2_top$TermCount))) +
  xlab("-log10(P.value)") +
  ylab("Function") +
  theme_minimal()

# Save
while (!is.null(dev.list()))  dev.off()
pdf("./results/Functional analysis TAM2/Bar_functions.pdf", width = 4, height = 4, family="ArialMT")
plotbarF
while (!is.null(dev.list()))  dev.off()


## Network plot
functions.TAM2_top.enr <- enrichDF2enrichResult(enrichDF = functions.TAM2_top, keyColname = "Term", 
                                                geneColname = "miRNA", pvalueColname = "P.value", 
                                                pvalueCutoff = 0.05, descriptionColname = "Term")
plotnetF <- cnetplot(functions.TAM2_top.enr, 
                     showCategory = 10,
                     node_label = "category",
                     cex_category = 1,
                     cex_gene = 1,
                     cex_label_category = 1,
                     cex_label_gene = 0.7,
                     color_category = "#E5C494",
                     color_gene = "#B3B3B3"#,
                     #hilight.params = list(category = c("Hematopoiesis"), alpha_hilight = 1, alpha_no_hilight = 0.3)
)

# Save
while (!is.null(dev.list()))  dev.off()
pdf("./results/Functional analysis TAM2/Net_functions.pdf", width = 7, height = 7, family="ArialMT")
plotnetF
while (!is.null(dev.list()))  dev.off()


### Diseases

diseases.TAM2 <- Enrichment.TAM2_result[Enrichment.TAM2_result$Category=="Disease",]
diseases.TAM2 <- diseases.TAM2[order(diseases.TAM2$Term, decreasing = F),]
diseases.TAM2 <- diseases.TAM2[order(diseases.TAM2$P.value, decreasing = F),]
diseases.TAM2_top <- diseases.TAM2[1:10,]

## Barplot
diseases.TAM2_top$TermCount <- paste0(diseases.TAM2_top$Term, " (", diseases.TAM2_top$Count, ")")
diseases.TAM2_top$logpval <- -log10(diseases.TAM2_top$P.value)
diseases.TAM2_top$TermCount <- factor(diseases.TAM2_top$TermCount, levels = unique(diseases.TAM2_top$TermCount))
library(ggplot2)
plotbarF <- ggplot(diseases.TAM2_top, aes(x=logpval, y=TermCount)) +
  geom_bar(stat = "identity") +
  scale_y_discrete(limits = rev(levels(diseases.TAM2_top$TermCount))) +
  xlab("-log10(P.value)") +
  ylab("Disease") +
  theme_minimal()

# Save
while (!is.null(dev.list()))  dev.off()
pdf("./results/Functional analysis TAM2/Bar_diseases.pdf", width = 4.8, height = 4, family="ArialMT")
plotbarF
while (!is.null(dev.list()))  dev.off()


## Network plot
diseases.TAM2_top.enr <- enrichDF2enrichResult(enrichDF = diseases.TAM2_top, keyColname = "Term", 
                                               geneColname = "miRNA", pvalueColname = "P.value", 
                                               pvalueCutoff = 0.05, descriptionColname = "Term")
plotnetD <- cnetplot(diseases.TAM2_top.enr, 
                     showCategory = 10,
                     node_label = "category",
                     cex_category = 1,
                     cex_gene = 1,
                     cex_label_category = 1,
                     cex_label_gene = 0.7,
                     color_category = "#E5C494",
                     color_gene = "#B3B3B3"#,
                     #hilight.params = list(category = c("Leukemia, Myeloid, Acute"), alpha_hilight = 1, alpha_no_hilight = 0.3)
)

# Save
while (!is.null(dev.list()))  dev.off()
pdf("./results/Functional analysis TAM2/Net_diseases.pdf", width = 7, height = 7, family="ArialMT")
plotnetD
while (!is.null(dev.list()))  dev.off()

#----
