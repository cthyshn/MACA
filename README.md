# MACA
An R package for multi-ancestry colocalization approaches. This package facilitates the integration of results from multi-ancestry fine-mapping methods, SuSiEx and MsCAVIAR, with colocalization methods, coloc and eCAVIAR. Four multi-ancestry colocalization approaches are available through the cross combinations of these methods: coloc_SuSiEx, eCAVIAR_SuSiEx, coloc_MsCAVIAR, eMsCAVIAR. 

## INSTALLATION 
```r
install.packages("devtools") 
devtools::install_github("cthyshn/MACA")
```

## STEP 1: MULTI-ANCESTRY FINE-MAPPING
This package relies on the multi-ancestry fine-mapping results from SuSiEx and/or MsCAVIAR. Follow the links below for the installation and usage instructions for your multi-ancestry fine-mapping method of choice. 
- [SuSiEx](https://github.com/getian107/SuSiEx)
- [MsCAVIAR](https://github.com/nlapier2/MsCAVIAR) <br>

Provide multiple ancestry-specific summary statistics for trait 1 to your selected method. Repeat for trait 2.

## STEP 2: COLOCALIZATION ANALYSIS 
#### coloc_SuSiEx and eCAVIAR_SuSiEx
After performing multi-ancestry fine-mapping of each trait (separately) using SuSiEx, provide the `.snp` files as input to `coloc_susiex()` or `eCAVIAR_SuSiEx()` for the colocalization step.
```r
trait1 <- fread("trait1.SuSiEx.EUR.AFR.output.snp")
trait2 <- fread("trait2.SuSiEx.EUR.AFR.output.snp")
coloc_susiex(trait1, trait2)
ecaviar_susiex(trait1, trait2)
```
#### coloc_MsCAVIAR and eMsCAVIAR
After performing multi-ancestry fine-mapping of each trait (separately) using MsCAVIAR, provide the `post.txt` files as input to `coloc_mscaviar()` or `emscaviar()` for the colocalization step. 
```r
trait1 <- fread("trait1.MsCAVIAR.EUR.AFR_post.txt")
trait2 <- fread("trait2.MsCAVIAR.EUR.AFR_post.txt")
coloc_mscaviar(trait1, trait2)
emscaviar(trait1, trait2)
```








