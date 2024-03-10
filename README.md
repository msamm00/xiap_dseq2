# xiap_dseq2

<br><br><br><br>
<b>This repository holds scripts for DSeq2 analysis of NCBI RNA sequencing data set GSE15168</b><br><br>
<br><br><br>
&emsp;<span>&#8226;</span> &nbsp; [NCBI GSE15168](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151648) &emsp;&emsp;&emsp; RNA Seq human liver data source 
with and without ischemic reperfusion injury measured at 2 different timepoints pre and post transplant<br><br>

&emsp;<span>&#8226;</span> &nbsp; [GSE15168 Raw Counts file](https://github.com/msamm00/xiap_dseq/GSE151648_liver-iri-counts.txt.gz) &emsp;&emsp; Raw counts file downloaded from NCBI data source<br><br>


&emsp;<span>&#8226;</span> &nbsp; [GSE15168 Raw Counts input file](https://github.com/msamm00/xiap_dseq/GSE151648_liver-iri-counts_cleaned.txt)  &emsp;&emsp;Raw counts file downloaded from NCBI data source cleaned up to match col names in the metadata file<br><br>



&emsp;<span>&#8226;</span> &nbsp; [GSE15168 SRA run info input file](https://github.com/msamm00/xiap_dseq/SraRunTable.txt)  &emsp; &emsp;&emsp;SRARun metadata information that contains details on what each sample means<br><br>



&emsp;<span>&#8226;</span> &nbsp; [runDESeq2.R](https://github.com/msamm00/xiap_dseq/runDESeq2.R)  &emsp;&emsp;&emsp; DESeq2 pipeline for running differential sequence expression analysis on RNA seq raw counts data<br><br>



| &emsp;<span>&#8226;</span> &nbsp; [postDESeq2.R](https://github.com/msamm00/xiap_dseq/postDESeq2.R)  &emsp; DESeq2 pipeline for running differential sequence expression analysis on RNA seq post transplant data<br><br>
