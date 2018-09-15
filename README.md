vicbits
=======

This is my collection of tools for DNA methylation analysis.

### Install

    library(devtools)
    install_github('wvictor14/vicbits')

Probably the only function that is general enough to be useful is the
`lmmatrix` function:

    library(vicbits) 

    # load example data
    library(minfiData) 
    data('RGsetEx')

    # calculate pcs on the data
    betas <- getBeta(RGsetEx)
    pc_obj <- prcomp(t(na.omit(betas)), center = T, scale = T)
    rotated <- pc_obj$x

    #rsquared
    lmmatrix(dep = rotated,
             ind = pData(RGsetEx)[,c('Sample_Group', 'age', 'sex', 'status')])

    ##                     PC1         PC2        PC3        PC4        PC5
    ## Sample_Group 0.73667341 0.007127756 0.21870722 0.03600302 0.00148859
    ## age          0.09045417 0.395734983 0.14766760 0.06044302 0.30570023
    ## sex          0.00625962 0.797344395 0.01945697 0.12347220 0.05346682
    ## status       0.73667341 0.007127756 0.21870722 0.03600302 0.00148859
    ##                    PC6
    ## Sample_Group 0.2498592
    ## age          0.5534485
    ## sex          0.1104967
    ## status       0.2498592

    #pvalue
    lmmatrix(dep = rotated,
             ind = pData(RGsetEx)[,c('Sample_Group', 'age', 'sex', 'status')],
             metric = 'Pvalue')

    ##                     PC1        PC2       PC3       PC4       PC5
    ## Sample_Group 0.02869701 0.87366182 0.3496483 0.7187988 0.9421553
    ## age          0.56246832 0.18086138 0.4519594 0.6386526 0.2551586
    ## sex          0.88157099 0.01657878 0.7921246 0.4946140 0.6593381
    ## status       0.02869701 0.87366182 0.3496483 0.7187988 0.9421553
    ##                     PC6
    ## Sample_Group 0.31265844
    ## age          0.08995476
    ## sex          0.51974948
    ## status       0.31265844
