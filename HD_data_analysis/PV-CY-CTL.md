Detailed protocol for the workflow PathVisio/Cytoscape/CyTargetLinker (PV-CY-CTL) used in the study "Network-based approaches for the interpretation of transcriptomics data in rare diseases, an application to Huntington’s". The RNAseq data was first produced and published by Labadorf et al. (PMID: 26636579) and the raw and processed data is available in GEO (GSE64810).

Input: 
* Differentially expressed genes from GSE64810
Software:
* PathVisio (version 3.0) PMID: 25706687
* Cytoscape (version 3.8.2) PMID: 14597658
* Cytoscape apps CyTargetLinker (version 4.0.0) and WikiPathways (version 3.3.7)
Datasets:
* WikiPathways human pathway collection (version 10.10.2020) [PMID: 33211851]
* CyTargetLinker linksets “Homo sapiens (hsa) - curated collection” and “Homo sapiens (hsa) - Reactome collection” release version 2021-01-10

Steps:
1. Install PathVisio, the WikiPathways app and the required databases for identifier mapping following these guidelines [https://docs.google.com/document/d/14Qe54Eb-Yz3RfLLDnBV3tRNwX0PcJOX5OLUknSMINzU/edit]. Creation of a WikiPathways account is not required for this analysis.
2. Import the RNAseq data to PathVisio following this tutorial [https://pathvisio.org/tutorials/multi-omics-tutorial.html]
3. Calculate pathway statistics as follows:
* Data -> Statistics
* Set the criterion to ([LogFC] < -1 OR [LogFC] > 1) AND [P.Value] ≤ 0.05
* Specify the pathway directory to the folder with the WikiPathways human pathway collection
* Calculate
