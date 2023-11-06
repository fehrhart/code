#Monte Carlo simulation test "large" - this test checks if groups of pathways overlap with a pathway when replaced randomly by genes available in WP. The results are given in 0 (no overlap) and 1 (overlapping genes (one or more!) found).
#clean workspace
rm(list=ls())
setwd("C:/...")

library(dplyr)

#1. import data
#import gene-pathway list derived from SPARQL query
gene_pathway <- read.csv("gene_pathway.csv", header= TRUE, fill = TRUE)  
as.data.frame.matrix(gene_pathway)

countlist <- read.csv("numbers_per.csv", header= TRUE, fill = TRUE)  
as.data.frame.matrix(countlist)

#subset list of available genes 
Available_genes <- (countlist$HGNC)

#2. create original gene subsets per pathway
WP4905 <- data.frame(c("ACP6","AFDN","AMELX","BCL9","CHD1L","CTNNB1","F11R","FMO5","GJA1","GJA3","GJA5","GJA8","KIRREL1","LINC00624","NBP12","OCLN","PRKAA1","PRKAA2","PRKAB1","PRKAB2","PRKAG1","PRKAG2","PRKAG3","PYGO1","PYGO2","TJP1","TJP2","TJP3"))
colnames(WP4905) <-("gene")
WP4906 <- data.frame(c("ADAM10","AKT1","BRINP1","CASP7","CEP19","CEP350","DLG1","DLG1-AS1","DYNC2H1","DYNC2LI1","DYNLL1","DYNLL2","DYNLRB1","DYNLRB2","DYNLT1","DYNLT3","FBXO45","FBXW7","FGFR1OP","FNDC8","GRIA1","HAMP","HFE","HIF1A","JUN","LINC00885","LINC01063","MAD2L1BP","MCRS1","MELTF","MELTF-AS1","MIR4797","MYC","MYCBP2","NCBP1","NCBP2","NCBP2-AS1","NCBP2-AS2","NF2","NRROS","PAK2","PCYT1A","PIGM","PIGX","PIGZ","PIK3R3","PXN","RABL2B","RABL2B","RNF168","RNF8","SDHAP1","SENP5","SIRT1","SLC40A1","SLC51A","SLC51B","SMCO1","STAT5A","STAT5B","TCTEX1D2","TF","TFRC","TGFB1","TM4SF19","TM4SF19-AS1","UBE2N","UBXN7","WDR34","WDR53","WDR60","ZDHHC19","ZNF76"))
colnames(WP4906) <-("gene")
WP4932 <- data.frame(c("ABHD11","ACACA","ACACB","ATF4","ATP5A1","ATP5B","ATP5D","ATP5E","ATP5F1","ATP5G1","ATP5G2","ATP5G3","ATP5O","ATPAF1","ATPAF2","BAZ1B","BCL7B","BECN1","BRD4","BTK","BUD23","CDKN1C","CFL1","CHTF18","CLASP1","CLASP2","CLDN1","CLDN3","CLDN4","CLDN5","CLIP2","CLTC","CTNNB1","DDX21","DEK","DLD","DLST","DNAJC30","EIF2A","EIF2AK3","EIF4H","ELN","ERCC6","ERCC6","FAS","FBLN2","FBLN5","FBN1","FKBP6","FZD9","GAPDH","GRB2","GRIP1","GTF2I","GTF2IRD1","H2AFX","HDAC2","HDAC3","HDAC6","HMGA1","HOXC8","HSPA2","LAT2","LIMK1","MAPK3","METTL27","MIR4284","MIR590","MLXIPL","MT-ATP6","MT-ATP8","MVB12A","MYBBP1A","MYC","MYO1C","NRG1","NUP62","OGDH","PCNA","PKLR","PRKG1","RB1","RFC2","RFC5","SF3B1","SMARCA5","SNAP25","SQSTM1","STX1A","TBL2","TMEM270","TRIM50","TSG101","UBE2E1","UBE2E3","UBE2L6","UBIAD1","ULK1","USF1","VAMP2","VPS28","VPS37A","VPS37B","VPS37C","VPS37D","VPS9D1","WNT2"))
colnames(WP4932) <-("gene")
WP4940 <- data.frame(c("CYFIP1","FMR1","NIPA1","NIPA2","TUBGCP2","TUBGCP3","TUBGCP4","TUBGCP5","TUBGCP6"))
colnames(WP4940) <-("gene")
WP4942 <- data.frame(c("CCL5","CHRFAM7A","CHRNA7","CREBBP","FAN1","FANCD2","FYN","GOLGA8H","GOLGA8Q","GOLGA8Q","GOLGA8R","GOLGA8R","GPR75","GRM6","KAT2B","KLF13","MIR211","MTMR10","OTUD7A","SERPINH1","TRPM1","ULK4P2"))   
colnames(WP4942) <-("gene")
WP4950 <- data.frame(c("ATP2A1","ATP2A1-AS1","ATXN2L","C3","CD19","CD81","CD82","CR2","GRB2","IFITM1","IL4","INSR","JAK2","KDR","LAT","MIR4721","MPL","NFATC2","NFATC2IP","PLN","PRMT1","RAB4A","RAB5A","RABEP2","RABGEF1","SDCCAG8","SH2B1","SLN","SPNS1","TRAF1","TRAF2","TUFM","VAV2"))
colnames(WP4950) <-("gene")
WP4949 <- data.frame(c("ALDOA","ASPHD1","BPTF","C16orf54","C16orf92","CASP8","CCDC6","CCT2","CCT3","CCT4","CCT5","CCT6A","CCT6B","CCT7","CCT8","CDIPT","CDIPT-AS1","CORO1A","DOC2A","ESR1","EZR","FAM57B","GDPD3","HDAC3","HIRA","HIRIP3","IGBP1","IGFBP3","INO80E","KCTD13","KIF22","KMT2C","KMT2D","MAP2K3","MAP2K6","MAPK3","MAZ","MIR3680-1","MIR3680-1","MIR3680-2","MIR3680-2","MSN","MVP","NFKB1","NR3C1","PAGR1","PARP4","PAXIP1","PCNA","PPARG","PPP2CA","PPP2CB","PPP2R1A","PPP2R5D","PPP4C","PPP4R1","PPP4R2","PPP4R3A","PPP4R3B","PPP4R3C","PPP4R4","PRRT2","PTEN","QPRT","REL","SEZ6L2","SIAH1","SPN","TAOK2","TBX6","TCP1","TMEM219","TP53","TRAF2","TRAF6","UNC13A","UNC13B","YPEL3","ZG16"))
colnames(WP4949) <-("gene")
WP4657 <- data.frame(c("SEPTIN5","SEPTIN8","SEPTIN11","ACTA2","ACTC1","AIFM3","ALDH1A2","ALDH4A1","ARNTL","ARVCF","ASF1A","BCL2","C22orf39","CBX5","CCDC188","CDC42","CDC45","CDH15","CHRD","CLDN1","CLDN3","CLDN5","CLTCL1","COMT","CRKL","CUL3","CYP26A1","CYP26B1","CYP26C1","DEPDC5","DGCR2","DGCR6L","DGCR8","DGCR9","DRD2","DROSHA","EGFR","EMC10","ESS2","FGF10","FGF8","FGFR1","FGFR2","FOXA2","FOXC1","FOXC2","GBX2","GLUD1","GNB1L","GP1BA","GP1BB","GP5","GP9","GSC2","HAND2","HDAC3","HES1","HIRA","HIRIP3","HIST1H4A","HIST1H4B","HIST1H4C","HIST1H4D","HIST1H4E","HIST1H4F","HIST1H4H","HIST1H4I","HIST1H4J","HIST1H4K","HIST1H4L","HIST2H4A","HIST2H4B","HIST4H4","KLHL22","KPNB1","LINC00896","LINC01637","LINC01662","LRRC74B","LZTR1","MAG","MALT1","MED15","MIR1306","MIR4761","MIR6816","MRPL40","NCOR1","NKX2-5","NPRL2","NPRL3","OAT","P2RX6","PAK4","PAX3","PI4KA","PITX2","PLK1","PPP1CB","PRKN","PRODH","RAF1","RAN","RANBP1","RANGAP1","RBX1","RCC1","RELN","RORC","RTL10","RTN4","RTN4R","SCARF2","SERPIND1","SHH","SHOC2","SLC25A1","SLC2A4","SLC7A4","SNAP29","SNORA77B","SREBF1","SREBF2","SRF","TANGO2","TBX1","THAP7","TMEM191A","TNPO1","TP53","TRMT2A","TSKS","TSSK2","TUBA3FP","TXNRD2","UFD1","VWF","XPO1","ZDHHC8","ZNF74"))
colnames(WP4657) <-("gene")
WP2380 <- data.frame(c("ACACB","ADAM17","AKT1","ALPL","APC","BAD","BCL2L11","BDNF","BMP2","CAMK1","CAMK2A","CAMK4","CASP3","CDC42","CDH2","CDK5","CDK5R1","CDKL5","CFL1","CHUK","CNR1","CREB1","CRTC1","CSNK2A1","CTNNB1","CYFIP1","DLG1","DOCK3","DOK5","DPYSL2","EEF2","EGR1","EGR2","EIF2S1","EIF2S2","EIF4E","EIF4EBP1","ELK1","FOS","FOXO3","FRS2","FRS3","FYN","GABRB3","GRB2","GRIA1","GRIA2","GRIA3","GRIN1","GRIN2B","GRIP1","GSK3B","HRAS","IGF2BP1","IKBKB","IKBKG","IRS1","IRS2","JAK2","JUN","KCNA3","KCNN2","KIDINS220","KSR1","LINGO1","MAP2K1","MAP2K2","MAP2K5","MAP3K1","MAP3K2","MAPK1","MAPK10","MAPK14","MAPK3","MAPK7","MAPK8","MAPK9","MAPT","MARCKS","MEF2A","MEF2C","MTOR","NCAM1","NCF1","NCF2","NCK1","NCK2","NFATC4","NFKB1","NFKBIA","NGF","NGFR","NSF","NTF3","NTRK1","NTRK2","NTRK3","PDPK1","PIK3R1","PIK3R2","PLCG1","PPP2CA","PRKAA1","PRKAA2","PRKCD","PTK2B","PTPN11","PTPRF","RAB3A","RAC1","RACK1","RAF1","RANBP9","RAP1A","RASGRF1","RELA","RHOG","RPS6","RPS6KA1","RPS6KA3","RPS6KA5","RPS6KB1","SH2B1","SH2B2","SHC1","SHC2","SHC3","SHC4","SIRPA","SORT1","SPP1","SQSTM1","SRC","STAT1","STAT3","STAT5A","STAT5B","SYN1","TIAM1","TRAF6","TSC2","VAV2","VAV3","YBX1"))
colnames(WP2380) <-("gene")

#___________________________________
#Repetition begin
#2. random sampling of genes according to their frequency in 8 pathways according to their number of genes
#unique genes per pathway, capped max freq
#counttable genes of WP - counttable genes in CNV - (remove those with very low freq which are only in CNV pathways) - adjust total number of genes in the respective CNV pathway

i <- 1
resultList = list()
repeat  {
#2. random sampling of genes
count_WP4905 <- dim(WP4905)[1]
n_WP4905 <- data.frame(sample(Available_genes, size = count_WP4905, replace = FALSE, prob = NULL))
colnames(n_WP4905) <-("gene")

count_WP4906 <- dim(WP4906)[1]
n_WP4906 <- data.frame(sample(Available_genes, size = count_WP4906, replace = FALSE, prob = NULL))
colnames(n_WP4906) <-("gene")

count_WP4932 <- dim(WP4932)[1]
n_WP4932 <- data.frame(sample(Available_genes, size = count_WP4932, replace = FALSE, prob = NULL))
colnames(n_WP4932) <-("gene")

count_WP4940 <- c(dim(WP4940)[1])
n_WP4940 <- data.frame(sample(Available_genes, size = count_WP4940, replace = FALSE, prob = NULL))
colnames(n_WP4940) <-("gene")

count_WP4942 <- c(dim(WP4942)[1])
n_WP4942 <- data.frame(sample(Available_genes, size = count_WP4942, replace = FALSE, prob = NULL))
colnames(n_WP4942) <-("gene")

count_WP4950 <- c(dim(WP4950)[1])
n_WP4950 <- data.frame(sample(Available_genes, size = count_WP4950, replace = FALSE, prob = NULL))
colnames(n_WP4950) <-("gene")

count_WP4949 <- c(dim(WP4949)[1])
n_WP4949 <- data.frame(sample(Available_genes, size = count_WP4949, replace = FALSE, prob = NULL))
colnames(n_WP4949) <-("gene")

count_WP4657 <- c(dim(WP4657)[1])
n_WP4657 <- data.frame(sample(Available_genes, size = count_WP4657, replace = FALSE, prob = NULL))
colnames(n_WP4657) <-("gene")

#3. compare if WP.... and WP2380 overlap, if dimension(row) of overlap > 0 set overlap to 1, sum up numbers at the end
WP2380_WP4905 <- inner_join(WP2380,n_WP4905)
if (dim(WP2380_WP4905)[1] == 0){WP2380_WP4905 <- c(0)}else{WP2380_WP4905 <- c(1)}

WP2380_WP4906 <- inner_join(WP2380,n_WP4906)
if (dim(WP2380_WP4906)[1] == 0){WP2380_WP4906 <- c(0)}else{WP2380_WP4906 <- c(1)}

WP2380_WP4932 <- inner_join(WP2380,n_WP4932)
if (dim(WP2380_WP4932)[1] == 0){WP2380_WP4932 <- c(0)}else{WP2380_WP4932 <- c(1)}

WP2380_WP4940 <- inner_join(WP2380,n_WP4940)
if (dim(WP2380_WP4940)[1] == 0){WP2380_WP4940 <- c(0)}else{WP2380_WP4940 <- c(1)}

WP2380_WP4942 <- inner_join(WP2380,n_WP4942)
if (dim(WP2380_WP4942)[1] == 0){WP2380_WP4942 <- c(0)}else{WP2380_WP4942 <- c(1)}

WP2380_WP4950 <- inner_join(WP2380,n_WP4950)
if (dim(WP2380_WP4950)[1] == 0){WP2380_WP4950 <- c(0)}else{WP2380_WP4950 <- c(1)}

WP2380_WP4949 <- inner_join(WP2380,n_WP4949)
if (dim(WP2380_WP4949)[1] == 0){WP2380_WP4949 <- c(0)}else{WP2380_WP4949 <- c(1)}

WP2380_WP4657 <- inner_join(WP2380,n_WP4657)
if (dim(WP2380_WP4657)[1] == 0){WP2380_WP4657 <- c(0)}else{WP2380_WP4657 <- c(1)}

#4. record individual overlaps
result <- c(WP2380_WP4905 + WP2380_WP4906 + WP2380_WP4932 + WP2380_WP4940 + WP2380_WP4942 + WP2380_WP4950 + WP2380_WP4949 + WP2380_WP4657)

resultList[[length(resultList)+1]] = result
i <- i+1
if (i>100) {break}
}
#_______________________________________________________
#repetition end
#this is set to i=100 repetitions for testing

#Results are in resultList, as double list, this step does type conversion and plotting of the results as mean, sd and histogram plot
resultList <- as.numeric(resultList)
mean(resultList)
sd(resultList)
hist(resultList)
