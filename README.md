DeconvBenchmark
===============
This repository holds all the data and code to perform our proposed benchmarking tests over immune cell deconvolution methods and molecular signatures. 
Performing these evaluations brings a framework for a fair measurement and comparision of the prediction accuracy and cell-type detection performance of any deconvolution method-molecular signature pair that the user would want to test. 

![1328](https://user-images.githubusercontent.com/68829569/191877734-4a5b1ad9-bc4d-4d3b-8f15-941fa8f8f365.png)

How to install and run the examples
------------------
Download the [zip file](http:://oasjdj.org) and unzip in your favorite directory in you local machine.

Clone this repository to your local device: "MyFavoriteDirectory/DeconvBenchmark".  
For a guided tutorial through the tests proposed in our benchmarking algorithm, load "DeconvBenchmark.Rproj" with RStudio to set your working environment and then run the file ["Tutorial.R"](https://github.com/agxnv/DeconvBenchmark/blob/main/Tutorial.R). 

Datasets
--------
Our tests requiere certain simulated and semi-simulated samples data to be performed. For our work, we used [murine sc-RNASeq data](https://figshare.com/articles/dataset/MCA_DGE_Data/5435866?file=10756795) from the **Mouse Cell Atlas** (_Han, X. et al. Mapping the Mouse Cell Atlas by Microwell-Seq. Cell, 2018_) to generate semi-simulated pooled samples DGE, wich are agregatted scRNAseq experiments from an unique cell type, simulating a pure bulk sample.

* **Bias test**: The user can generate simulated bulk samples based on any molecular signature to test with the function [SimulatedMixtures](https://rdrr.io/github/elmerfer/MIXTURE/man/SimulatedMixtures.html) from the [MIXTURE](https://github.com/elmerfer/MIXTURE) package.
* **Null test**: For semi-simulated pooled samples DGE containing parenchymatous and stromal cells from murine tissues, download the file ["simulated_nonimmune.RData"](https://github.com/agxnv/DeconvBenchmark/blob/main/Data/simulated_nonimmune.RData).
* **True positive cases test**: For semi-simulated pooled samples DGE containing just one immune cell-type cells per case from murine tissues, download the file ["simulated_immune.RData"](https://github.com/agxnv/DeconvBenchmark/blob/main/Data/simulated_immune.RData).

All the semi-simulated pooled samples were transformed from raw counts to TPM. For calculating the length of each gene, we took the reference genome used to annotate the expression matrices ([Mus musculus GRCm 38.88](https://ftp.ensembl.org/pub/release-88/gtf/mus_musculus/Mus_musculus.GRCm38.88.gtf.gz)), aggregated and averaged the exon’s length for each version of that gene available. See the file ["tpm_calculator.R"](https://github.com/agxnv/DeconvBenchmark/blob/main/Data/tpm_calculator.R).

Molecular signatures
--------------------
We tested the [murine tissue-specific molecular signatures](https://github.com/wuaipinglab/ImmuCC/tree/master/tissue_immucc/SignatureMatrix) developed by **Chen et al.** (_Ziyi Chen, Chengyang Ji, Qin Shen, Wei Liu, F Xiao-Feng Qin, Aiping Wu, Tissue-specific deconvolution of immune cell composition by integrating bulk and single-cell transcriptomes, Bioinformatics, 2020_). 

Deconvolution methods
---------------------
We tested six reference-based methods made to deconvolve bulk RNA-Seq samples. 

| Method | Version | Reference |
|--------|---------|-----------|
| **CIBERSORTx** | Web ([cibersortx.stanford.edu](https://cibersortx.stanford.edu)) | Newman, A.M., Steen, C.B., Liu, C.L. et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. Nat Biotechnol, 2019 |
| **DCQ** | R package: [ADAPTS v1.0.21](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0224693) | Altboum, Z. et al. Digital cell quantification identifies global immune cell dynamics during influenza infection. Mol. Syst. Biol. 10, 720, 2014 |
| **DeconRNASeq** |  R package: [ADAPTS v1.0.21](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0224693) | Gong, T. & Szustakowski, J. D. DeconRNASeq: a statistical framework for deconvolution of heterogeneous tissue samples based on mRNA-Seq data. Bioinforma. Oxf. Engl. 29, 1083–1085, 2013 |
| **EPIC** | R package: [EPIC v1.1.5](https://github.com/GfellerLab/EPIC) |Racle, J., Jonge, K. de, Baumgaertner, P., Speiser, D.E., and Gfeller, D. Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. eLife, 6, e26476, 2017 |
| **MIXTURE** | R package: [MIXTURE v0.0.1](https://github.com/elmerfer/MIXTURE) | Fernández EA, Mahmoud YD, Veigas F et al. Unveiling the immune infiltrate modulation in cancer and response to immunotherapy by MIXTURE-an enhanced deconvolution method. Brief Bioinform, 20;22(4), 2021 |
| **quanTIseq** | R package: [immunedeconv v2.0.4](https://link.springer.com/protocol/10.1007/978-1-0716-0327-7_16) | Finotello F, Mayer C, Plattner C, Laschober G, Rieder D, Hackl H, Krogsdam A, Loncova Z, Posch W, Wilflingseder D, Sopper S, Ijsselsteijn M, Brouwer TP, Johnson D, Xu Y, Wang Y, Sanders ME, Estrada MV, Ericsson-Gonzalez P, Charoentong P, Balko J, de Miranda NFDCC, Trajanoski Z. Molecular and pharmacological modulators of the tumor immune contexture revealed by deconvolution of RNA-seq data. Genome Medicine, 2019 |

Authors
-------
* **Elmer A. Fernández** - *Idea* - [Profile](https://www.researchgate.net/profile/Elmer_Fernandez) - [CIDIE]- [CONICET]
* **Agustín Nava** (PhD candidate) - [agustin.nava@huesped.org.ar](mailto:agustin.nava@huesped.org.ar)
