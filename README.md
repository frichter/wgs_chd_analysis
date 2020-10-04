## Genomic analyses implicate noncoding de novo variants in congenital heart disease

**Contributors:** Felix Richter (felix.richter@icahn.mssm.edu), Sarah Morton (sarah.morton@childrens.harvard.edu), Alex Kitaygorodsky (ak3792@cumc.columbia.edu), Kathleen Chen (kchen@flatironinstitute.org)

**DNV identification:** source code is available at https://github.com/frichter/dnv_pipeline (GATK preprocessing and FreeBayes) and https://github.com/ShenLab/igv-classifier (for IGV-driven neural network methods).

**HeartENN:** Please refer to the documentation in our new [HeartENN-models](https://github.com/FunctionLab/HeartENN-models) repository. Note that HeartENN-models is currently a work-in-progress; we are still building out the scripts and documentation to make HeartENN easier to run. In the meantime, please reach out to kc31@princeton.edu for more information.  

**Burden tests:** Scripts used for key steps in WGS de novo variant burden tests are included in this repository.
- RBP_burden_analysis.R: burden of variants in post-transcriptional regulatory regions
- multihit_enhancer_identification.R: over-representation of genes with multiple noncoding DNVs in prioritized enhancers
- heartenn_burden.R: burden of high-scoring HeartENN variants in cases vs controls

**Implementation:** Please note that this repository is intended as a means of providing code review, and the encolsed scripts may require additional local dependencies. More extensive custom scripts and pipelines are available from the authors on request.

**Relevant data:** Supplemental data can be found in supplementary_tables_chd_wgs.xlsx

