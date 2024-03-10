# xiap_dseq2

<br><br><br><br>
<b>This repository holds scripts for DSeq2 analysis of NCBI RNA sequencing data set GSE15168</b><br><br>
<br><br><br>

| Name | Description |
| --- | --- |
|[NCBI GSE15168](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151648) | RNA Seq human liver data source with and without ischemic reperfusion injury measured at 2 different timepoints pre and post transplant<br><br> |
| [GSE15168 Raw Counts file](https://github.com/msamm00/xiap_dseq/GSE151648_liver-iri-counts.txt.gz) | Raw counts file downloaded from NCBI data source |
| [GSE15168 Raw Counts input file](https://github.com/msamm00/xiap_dseq/GSE151648_liver-iri-counts_cleaned.txt)  | Raw counts file downloaded from NCBI data source cleaned up to match col names in the metadata file |
| [GSE15168 SRA run info input file](https://github.com/msamm00/xiap_dseq/SraRunTable.txt)  | SRARun metadata information that contains details on what each sample means |
| [runDESeq2.R](https://github.com/msamm00/xiap_dseq/runDESeq2.R)  | DESeq2 pipeline for running differential sequence expression analysis on RNA seq raw counts data |
|[postDESeq2.R](https://github.com/msamm00/xiap_dseq/postDESeq2.R) | DESeq2 pipeline for running differential sequence expression analysis on RNA seq post transplant data |
