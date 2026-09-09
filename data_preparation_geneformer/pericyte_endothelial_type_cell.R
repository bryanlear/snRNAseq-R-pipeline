library(Seurat)
library(Matrix)

object <- readRDS("results/GSE174367/06_ambient_rna/Sample-22.rds")

metadata <- object[[]]
barcodes <- rownames(metadata)[which(metadata$Cell.Type == "PER.END")]

before <- LayerData(object, assay = "RNA", layer = "counts")[, barcodes, drop = FALSE]
after <- LayerData(object, assay = "decontX", layer = "counts")[, barcodes, drop = FALSE]
per_end_metadata <- metadata[barcodes, , drop = FALSE]

counts_before <- colSums(before)
counts_after <- colSums(after)

print(data.frame(barcode = colnames(before), counts_before = counts_before, 
counts_after = counts_after, percent_removed = 100 * (counts_before - counts_after) / counts_before,row.names = NULL))

#### SAMPLE 19
#               barcode counts_before counts_after percent_removed
# 1  AAGAACATCTTAGCAG-1          4803         4772      0.64542994
# 2  AAGCATCCACTGTCCT-1          9151         9151      0.00000000
# 3  ACCATTTAGGTAAAGG-1          4479         4479      0.00000000
# 4  AGTCACACATCGCCTT-1          6450         6447      0.04651163
# 5  AGTTCGACATAGAATG-1          6805         6799      0.08817046
# 6  ATCCGTCCAATCTGCA-1          9289         9289      0.00000000
# 7  ATTCGTTAGGAGCTGT-1          1784         1764      1.12107623
# 8  ATTGTTCGTGCAACGA-1          6062         6062      0.00000000
# 9  CACGGGTCACTTGAGT-1          2951         2661      9.82717723
# 10 CTCCTCCCAGGCACAA-1          2555         2509      1.80039139
# 11 GAAGAATAGAGGGTAA-1         11990        11990      0.00000000
# 12 GGTGTCGAGAGCCGTA-1          5325         5322      0.05633803
# 13 GTACAACTCTTAATCC-1          8596         8596      0.00000000
# 14 GTGTTCCCACTAGTAC-1          4971         4970      0.02011668
# 15 GTTGTGATCAAAGGTA-1          6569         6518      0.77637388
# 16 TCTCAGCTCCTTTGAT-1          9383         9383      0.00000000
# 17 TGAGTCAGTATTCCGA-1          5778         5775      0.05192108
# 18 TGTGGCGAGGTTACAA-1          3224         3201      0.71339950
# 19 TGTTTGTTCGACCATA-1          1821         1798      1.26304228
# 20 TTCAATCTCCCTCTCC-1          3613         3612      0.02767783


#### SAMPLE 22

#              barcode counts_before counts_after percent_removed
# 1 GAGCTGCTCTTCGTGC-4          2167         1448        33.17951
# 2 GGGAGTATCTCCGCAT-4          1884          734        61.04034
# 3 GTGATGTAGGGACACT-4          1649          603        63.43238
# 4 TCCACCATCACTCTTA-4          1465          491        66.48464
# 5 TCGCTTGTCCGTTGAA-4          2773         1119        59.64659
# 6 TGCTCGTGTAAGGCCA-4          5546         2930        47.16913