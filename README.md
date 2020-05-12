## Genomic analyses implicate noncoding de novo variants in congenital heart disease

**Contributors:** Felix Richter (felix.richter@icahn.mssm.edu), Sarah Morton (sarah.morton@childrens.harvard.edu), Alex Kitaygorodsky (ak3792@cumc.columbia.edu), Kathleen Chen (kchen@flatironinstitute.org)

**DNV identification:** source code is available at https://github.com/frichter/dnv_pipeline (GATK preprocessing and FreeBayes) and https://github.com/ShenLab/igv-classifier (for IGV-driven neural network methods).

**HeartENN:** Model weights used for the variant effect prediction are available at [HeartENN_models.tar.gz](https://www.dropbox.com/s/rerhx0sdc6y8r6v/HeartENN_models.tar.gz?dl=0) (Dropbox link, 865 Mb). These model weights, trained on 184 cardiac epigenomic annotations, can be used in conjunction with [Selene](https://github.com/FunctionLab/selene/archive/0.4.8.tar.gz) to predict variant effects.

**Burden tests:** Scripts used for key steps in WGS de novo variant burden tests are included in this repository.
- RBP_burden_analysis.R: burden of variants in post-transcriptional regulatory regions
- multihit_enhancer_identification.R: over-representation of genes with multiple noncoding DNVs in prioritized enhancers
- heartenn_burden.R: burden of high-scoring HeartENN variants in cases vs controls

**Implementation:** Please note that this repository is intended as a means of providing code review, and the encolsed scripts may require additional local dependencies. More extensive custom scripts and pipelines are available from the authors on request.

**Relevant data:** Supplemental data can be found in supplementary_tables_chd_wgs.xlsx

