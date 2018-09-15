# vicbits

This is my collection of tools for DNA methylation analysis.

### Install

```r
library(devtools)
install_github('wvictor14/vicbits')
```

Probably the only function that is general enough to be useful is the `lmmatrix` function:

```r
library(vicbits)
library(minfiData)

data('RGsetEx')

betas <- getBeta(RGsetEx)
pc_obj <- prcomp(t(na.omit(betas)), center = T, scale = T)
rotated <- pc_obj$x

lmmatrix(dep = rotated,
         ind = pData(RGsetEx)[,c('Sample_Group', 'age', 'sex', 'status')])
         
# p value

lmmatrix(dep = rotated,
         ind = pData(RGsetEx)[,c('Sample_Group', 'age', 'sex', 'status')],
         metric = 'Pvalue')
```