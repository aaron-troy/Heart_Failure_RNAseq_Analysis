# Heart_Failure_RNAseq_Analysis
Analysis of gene expression profiles obtained via RNA sequencing of cardiac biopsies obtained from heart failure patients and healthy controls. Data is associated with the publication "Myocardial Gene Expression Signatures in Human Heart Failure With Preserved Ejection Fraction" by Hahn et al, Circulation 2020. 

## Differentially Expressed Gene (DEG) Extraction

There is inevitably a lot of variation introduced to the data by differences in tissue collection, handling, processing, etc. sample to sample. Variabillity in the presence in contaminating blood is one of many components here. To suppress this, I removed genes with a correlation > 0.5 with the β-globin component of hemoglobin. A real valued paitent variables (e.g. BMI) were converted to Z scores prior to any modelling. 

To complete count normalization and indentification of differentially expressed genes (DEGs), I leveraged the R package DESeq2,  

#### Normalization and Variance Stabillization
In the RNA sequencing workflow there are a number of distinctions between samples that are often present in raw counts data. These include: 
- Sample to sample differences in library size and hence total number of reads. 
- Differences in the library composition. For example, a gene highly expressed in one sample may be
 knocked out in another. Reads attributed to this gene in the WT sample will be distributed to others in the KO, which will appear as if many gene are DE, but this is clearly not the case.

There are several normalizing approaches that can be applied to account for these differences. Here I applied the DESeq2 implementation of median-of-ratios normalization. 

This was followed by variance stabillization, a necessary step to deal with the high degree of [heteroscedasticity](https://en.wikipedia.org/wiki/Homoscedasticity_and_heteroscedasticity) typically present in sequencing data. This is a prerequsite for most statistical analyses you'd want to apply (regression, ANOVA, etc.) - anything that expects variance in a gene to be independant of mean expression. 

#### Effect of Covariates
Something notable within this initial exploratory analysis was the clear effect of sex in the counts. Given the interest in HFpEF for this project, I adjusted for this variable prior. 
![alt text](image.png)
After the adjustment the effect of sex if absent, and the presence of three clusters (control, HFpEF, HFrEF) is clear. 
![alt text](image-1.png)
![alt text](image-2.png)

#### Differentially Expressed Genes
![alt text](image-3.png)

#### Enriched GO Biological Processes
![alt text](image-4.png)
![alt text](image-5.png)

## Transcription Factor Inference using VIPER 
![alt text](image-6.png)
![alt text](image-7.png)
## Exploring Protein-Protein Interactions through Prize Collecting Steiner Forests 

