# miRNA-predictive-signature-survival

Pipeline for creating and evaluating a survival predictive miRNA signature using miRNA expression data.

## Study

A comprehensive review of all previous studies that have analyzed the expression of microRNAs as predictors of survival in pediatric acute myeloid leukemia (AML) was conducted. Subsequently, we utilized miRNA expression data from a publicly available cohort of pediatric AML patients to develop a novel predictive model and compared the predictive accuracy of our signature with all previous models. Herein, we propose a novel score based on a 37-microRNA signature that is associated with AML and is able to predict survival with greater accuracy than previous microRNA-based methods.

## Data

We used miRNA expression and clinical data from pediatric AML patients that is freely available through the Therapeutically Applicable Research to Generate Effective Treatments (https://ocg.cancer.gov/programs/target) initiative, at the Genomic Data Commons repository (https://portal.gdc.cancer.gov/projects). For signature quality control we used an external gene expression dataset, publicly available at the Gene Expression Omnibus repository (GSE97135).

## Workflow

To follow the pipeline, run the R scripts in `/src` in the order indicated.
As randomization was employed in certain sections of the original analysis, some results are already available in the `/results` directory for users to replicate the original article's findings.

## Reference

Article under review

