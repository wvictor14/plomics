---
title: "plomics"
output:
  html_document:
    df_print: kable
    keep_md: yes
    self_contained: yes
    theme: spacelab
    toc: yes
    toc_depth: 3
    toc_float:
      collapsed: no
editor_options:
  chunk_output_type: console
---

A collection of functions for placental DNA methylation analysis.

- [Install](#install)
  - [lmmatrix](#lmmatrix)
  - [pairTest](#pairtest)
  - [findsentrix](#findsentrix)


## Install


```r
remotes::install_github('wvictor14/plomics')
```



## Functions

### lmmatrix

Computes pairwise linear models between several variables.


```r
library(minfiData)
library(plomics)

# load example data
data(RGsetEx)

# calculate pcs on the data
betas <- getBeta(RGsetEx)
pc_obj <- prcomp(t(na.omit(betas)), center = T, scale = T)

# get pc scores for each sample
rotated <- pc_obj$x

#rsquared
rsq <- lmmatrix(dep = rotated,
                ind = as.data.frame(pData(RGsetEx)[,c('Sample_Group', 'age', 'sex', 'status')]))

#pvalue
pva <- lmmatrix(dep = rotated,
                ind = as.data.frame(pData(RGsetEx)[,c('Sample_Group', 'age', 'sex', 'status')]),
                metric = 'Pvalue')
##### plot
# reshape first
rsq_plot <- rsq %>% as.data.frame() %>% 
  
  # add dep variables
  mutate(dep = rownames(rsq)) %>%
  
  # reshape
  gather(PC, rsquared, -dep)

pva_plot <- pva %>% as.data.frame() %>% 
  
  # add dep variables
  mutate(dep = rownames(rsq)) %>%
  
  # reshape
  gather(PC, pval, -dep) %>%
  
  # pvalue categories
  mutate(pval_cat = case_when(
    pval > 0.05  ~ '> 0.05',
    pval < 0.05 & pval > 0.01 ~ '< 0.05',
    pval < 0.01 & pval > 0.001 ~ '< 0.01',
    pval < 0.001 ~ '< 0.001'
  ))
  
ggplot(rsq_plot, aes(x = PC, y = dep, fill = rsquared)) +
  geom_tile() + theme_bw() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colours=c("white", "#ffffcc", "#41b6c4", "#2c7fb8", "#253494"), 
                       breaks = c(0,0.5,1), limits = c(0,1), 
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) 
```

![](README_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
ggplot(pva_plot, aes(x = PC, y = dep, fill = pval_cat)) +
  geom_tile() + theme_bw() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = c('> 0.05' = 'white', '< 0.05' = '#fee8c8', 
                               '< 0.01' = '#fdbb84', '< 0.001' = '#e34a33'))
```

![](README_files/figure-html/unnamed-chunk-3-2.png)<!-- -->

### pairTest

To test if covariates are confounding each other, we need to pairwise tests of independence between
covariates. If at least one covariate is numeric, we can use linear regression. Otherwise if both
covariates are categorical, then a chi squared test must be used.


```r
variables <- data.frame(
  row = c(1, 1, 2, 2, 3, 3),
  column = c(1, 1, 1, 2, 2, 2),
  Sex = c('m', 'm', 'm', 'f', 'f', 'f'),
  age = c(18, 19, 18, 27, 30, 16),
  ethnicity = c('AF', 'AF', 'AF', 'EU', 'EU', 'AS')
)

tests <- pairtest(variables)
# make categories
tests <- tests %>% 
  mutate(pval_cat = if_else(p.value < 0.001, '< 0.001',
                            if_else(p.value < 0.01, '< 0.01',
                                    if_else(p.value < 0.05, '< 0.05', '<1'))))
tests
```

<div class="kable-table">

X1   Row      Column              Fstat   df   p.value   Chi.Square  pval_cat 
---  -------  ----------  -------------  ---  --------  -----------  ---------
1    row      column       8.000000e+00    4     0.047           NA  < 0.05   
2    row      Sex          8.000000e+00    4     0.047           NA  < 0.05   
3    row      age          5.660000e-01    4     0.494           NA  <1       
4    row      ethnicity    3.643000e+00    3     0.104           NA  <1       
5    column   Sex          1.030218e+31    4     0.000           NA  < 0.001  
6    column   age          1.976000e+00    4     0.233           NA  <1       
7    column   ethnicity    7.605904e+30    3     0.000           NA  < 0.001  
8    Sex      age          1.976000e+00    4     0.233           NA  <1       
9    Sex      ethnicity              NA    2     0.050            6  <1       
10   age      ethnicity    4.591900e+01    3     0.221           NA  <1       

</div>

```r
# plot heatmap of associations
ggplot(tests, aes(x=Row, y = Column, fill = pval_cat)) +
  geom_tile(col = 'grey') + theme(panel.background = element_blank()) + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = c('> 0.05' = 'white', '< 0.05' = '#fee8c8', 
                               '< 0.01' = '#fdbb84', '< 0.001' = '#e34a33'))
```

![](README_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

### findsentrix

`findsentrix` takes a vector of sentrix IDs (chip identifiers) and searches a directory for IDAT 
files that match. The returned data frame contains two columns: (1) the sentrix ID (2) unique file 
paths for each idat that matches the sentrix ID. 

Here is an example using the robinson lab master sample sheet:


```r
# read in master sample sheet
ss <- readxl::read_xlsx('Z:/ROBLAB6 Infinium450k John/Master_Sample_Sheet.xlsx')

## specify idat directory
idat_dir <- 'Z:/ROBLAB6 Infinium450k John/EPIC Raw data/'

ss <- ss %>% 

# Take the first 6 EPIC samples
  dplyr::arrange(desc(Platform)) %>% 
  dplyr::slice(1:6) %>% 
  dplyr::select(Sample_Name, Sentrix_ID, Sentrix_Position) %>%
 
# create sentrix column
  dplyr::mutate(Sentrix = paste0(Sentrix_ID, '_', Sentrix_Position))
```

We created a sentrix ID by taking the chip serial number (confusingly named as "sentrix_ID") and 
pasting this to the position identifier ("Sentrix_Position"):

**Sentrix:** 200889820007_R01C01

**Sentrix Chip Number:** 200889820007

**Sentrix Position on chip:** R01C01


Now we can use `findsentrix` to find the filepaths for idats that match each sentrix identifier:


```r
idatfiles <- findsentrix(sentrix = ss$Sentrix, directory = idat_dir)
idatfiles
```

<div class="kable-table">

Sentrix               Basename                                                                                   
--------------------  -------------------------------------------------------------------------------------------
200889820007_R01C01   Z:/ROBLAB6 Infinium450k John/EPIC Raw data//Batch7_rescan/200889820007/200889820007_R01C01 
200889820007_R02C01   Z:/ROBLAB6 Infinium450k John/EPIC Raw data//Batch7_rescan/200889820007/200889820007_R02C01 
200889820007_R03C01   Z:/ROBLAB6 Infinium450k John/EPIC Raw data//Batch7_rescan/200889820007/200889820007_R03C01 
200889820007_R04C01   Z:/ROBLAB6 Infinium450k John/EPIC Raw data//Batch7_rescan/200889820007/200889820007_R04C01 
200889820007_R05C01   Z:/ROBLAB6 Infinium450k John/EPIC Raw data//Batch7_rescan/200889820007/200889820007_R05C01 
200889820007_R06C01   Z:/ROBLAB6 Infinium450k John/EPIC Raw data//Batch7_rescan/200889820007/200889820007_R06C01 

</div>

```r
# join all matches, retaining unmatched and multiple matched IDs
ss <- ss %>%
  dplyr::full_join(idatfiles, by = 'Sentrix')
```

Finally, we can load idats using this dataframe:


```r
## Now you can load in these samples with minfi::read.metharray.exp
rgset <- minfi::read.metharray.exp(targets = as.data.frame(ss), verbose = T)
```

```
## [read.metharray] Reading 200889820007_R01C01_Grn.idat
```

```
## [read.metharray] Reading 200889820007_R02C01_Grn.idat
```

```
## [read.metharray] Reading 200889820007_R03C01_Grn.idat
```

```
## [read.metharray] Reading 200889820007_R04C01_Grn.idat
```

```
## [read.metharray] Reading 200889820007_R05C01_Grn.idat
```

```
## [read.metharray] Reading 200889820007_R06C01_Grn.idat
```

```
## [read.metharray] Reading 200889820007_R01C01_Red.idat
```

```
## [read.metharray] Reading 200889820007_R02C01_Red.idat
```

```
## [read.metharray] Reading 200889820007_R03C01_Red.idat
```

```
## [read.metharray] Reading 200889820007_R04C01_Red.idat
```

```
## [read.metharray] Reading 200889820007_R05C01_Red.idat
```

```
## [read.metharray] Reading 200889820007_R06C01_Red.idat
```

```
## [read.metharray] Read idat files in 2.5 seconds
```

```
## [read.metharray] Creating data matrices ... done in 10.8 seconds
## [read.metharray] Instantiating final object ... done in 0.1 seconds
```

```r
rgset
```

```
## class: RGChannelSet 
## dim: 1052641 6 
## metadata(0):
## assays(2): Green Red
## rownames(1052641): 1600101 1600111 ... 99810990 99810992
## rowData names(0):
## colnames(6): 200889820007_R01C01 200889820007_R02C01 ...
##   200889820007_R05C01 200889820007_R06C01
## colData names(6): Sample_Name Sentrix_ID ... Basename filenames
## Annotation
##   array: IlluminaHumanMethylationEPIC
##   annotation: ilm10b4.hg19
```

