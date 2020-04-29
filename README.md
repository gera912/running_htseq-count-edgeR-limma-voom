# RNA-Seq differential expression analysis
The objectives for this assignment include counting aligned RNA-Seq reads across gene models for several libraries, normalization of RNA-Seq count data, and differential gene expression analysis using edgeR and limma with voom.

The data: You will be analyzing data from an experiment in which we exposed larval stickleback fish from two populations (Boot Lake “Bt” and Rabbit Slough “RS”) to one of two treatments: a conventional microbiota (“CV”) or germ-free (“GF”). RNA was isolated from the intestinal tracts of the experimental fish and used to make individual TruSeq mRNA-seq libraries. We are interested in genes differentially expressed between fish of different populations and genes differentially expressed between fish exposed or unexposed to microbes. We also might want to identify genes that respond to the treatment differently depending on a fish’s population (a.k.a. an interaction between population and treatment), but we’ll leave that for the fall.

All files for this assignment may be found in: /projects/bgmp/shared/Bi623/assign7/
Useful resources: Never neglect the usefulness of the help() function in R! edgeR user guide: 

http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf limma user guide: 

http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
# 1. HTSeq-count
Write an sbatch script (submitted as a job on Talapas) to count the number of reads mapping uniquely to each gene in the stickleback reference, for all SAM files.

These SAM files (which begin with “Bt_”) are alignments for 6 stickleback gut RNA-Seq libraries.

You can use htseq-count https://htseq.readthedocs.io/en/release_0.11.1/count.html?highlight=htseq-count#counting-reads-in-features-with-htseq-count to accomplish this.

Load the Talapas module HTSeq to run htseq-count as follows: ml purge ml slurm easybuild intel/2017a HTSeq/0.9.1-Python-3.6.1

Review the above link and the slide from lecture to run htseq-count, specifying that you do not have stranded data, that you want to restrict counting to exon alignments, and try the “union” criteria for counting reads. You can name your output file Gacu_readcounts.txt (10 points)

Note the column totals calculated by htseq-count, at the bottom of your file. (Naturally, you would normally exclude these in R before any kind of differential expression analysis.) What information do these totals contain? (5 points)

# 2. Evaluating the effect of normalization on count data

In R, read in the file /projects/bgmp/shared/Bi623/assign7/Gacu_gut_counts.tsv (name the df “Counts”). This is NOT your output from part 1. These are the gene-wise counts (number of reads aligning to gene models in the annotation) for a number of stickleback individuals. When you read in the file, set row.names=1, and stringsAsFactors=FALSE.

Using the unnormalized values in Counts, plot library 1A_02 values vs. library 1A_03 values. Draw a line through the origin (intercept=0) with slope=1, which would be the expectation for equal counts between the two libraries.

Hint: use plot() and abline() functions. What do you notice about the relationship between counts for these libraries? (5 points)

Next, store the raw counts as a “DGEList” object (remember, you’ll have to load the edgeR package, which you should already have installed.)
```
dge <- DGEList(counts=Counts)
```
Calculate TMM normalization factors and add to dge

```
dge <- calcNormFactors(dge)
```

Apply the TMM normalization factors and adjust for library size by converting to copies per million (cpm)
```
TMMCPM <- as.data.frame(cpm(dge,
                        normalized.lib.sizes=TRUE))
```
Now plot the normalized expression values for library 1A_02 vs those for library 1A_03 as above, including the slope=1 line. How has the relationship changed, and what does it mean regarding the effectiveness of normalization? (5 points)

# 3. Testing and visualizing effects of Population and Treatment on differential expression using the negative binomial generalized linear model approach in edgeR (20 points total)

First, we will define a new DGEList object to work with
```
dge.er <- DGEList(counts=Counts)
```
Now get rid of genes expressed at low levels in many of the samples, which are not reliably tested. We only keep genes with at least one read per million (cpm) in at least 8 samples.
```
keep.dge.er <- rowSums(cpm(dge.er)>=1) >= 8
dge.er <- dge.er[keep.dge.er,]
```
How many genes are we dealing with now? (5 points)

Now calculate TMM normalization factors
```
dge.er <- calcNormFactors(dge.er)
```

Now we set up the model using our experimental design, but first read in the Gacu_gut_metadata.tsv file using read.delim(), and name it “Targets.er”

One easy way to set up the model is to combine our Population and Treatment variables into one. This will allow us flexibility when we set up “contrasts” to test. Don’t worry about contrasts too much now. You’ll cover them in the biostatistics course. The general idea here is that there are 4 Population-Treatment factor level combinations (Bt-CV, Bt-GF, RS-CV, RS-GF). You can use contrasts to test for differences in gene expression between subsets of these groupings

```
PopTrt <- paste(Targets.er$Population,
                Targets.er$Treatment, sep=".")
PopTrt
```
Now we use the model.matrix() function to specify the experimental design. This shorthand way of doing it will set up a generalized linear model that includes the “intercept” (0) and the effects of the 4 factor level combinations.
```
design.er <- model.matrix(~0 + PopTrt)
```
Notice that we have 4 groups: Bt.CV, Bt.GF, RS.CV, and RS.GF
```
design.er
colnames(design.er)
```
We can now use these groups in contrasts to test hypotheses of interest

Before we conduct the hypothesis tests, we need to estimate the dispersion parameter for each gene. We serially estimate common, trended, and tagwise dispersions for design.er. Remember, this is necessary for modeling the dispersion parameter of the negative binomial distribution.
```
dge.er <- estimateGLMCommonDisp(dge.er, design.er)
dge.er <- estimateGLMTrendedDisp(dge.er, design.er)
dge.er <- estimateGLMTagwiseDisp(dge.er, design.er)
```
Now we fit the full generalized linear model, incorporating the dispersion estimates.
```
fit.er <- glmFit(dge.er, design.er)
```
Now we perform likelihood ratio tests for factor level contrasts of interest. These effectively test whether the likelihood of the data given a model with a particular effect term (e.g. Population) is higher than the likelihood of the data given the model without that effect term.

Population contrast To test for the Population effect we contrast the average of Bt.CV and Bt.GF vs. the average of RS.CV and RS.GF
```
lrtPop <- glmLRT(fit.er, contrast=c(-.5,-.5,.5,.5))
```
To look at the genes with the lowest FDR-corrected p-values:
```
topTags(lrtPop)
```
And to write all of the test results to file
```
write.table(topTags(lrtPop, n=16877), 'edgeR_PopLRT.txt',
            sep='\t')
```
Summarize how many genes are differentially expressed by Population

```
de.er <- decideTestsDGE(lrtPop, p=0.05)
summary(de.er)
```
Visualize the differential expression patterns with a “smear plot”
```
detags.er <- rownames(fit.er)[as.logical(de.er)]
plotSmear(lrtPop, de.tags=detags.er)
abline(h = c(-2, 2), col = "blue")
```
What do the red dots represent, and what boundaries do the blue lines mark?

Treatment Contrast To test for the Treatment (microbiota) effect we contrast the average of Bt.CV and RS.CV vs. the average of Bt.GF and RS.GF
```
lrtTrt <- glmLRT(fit.er, contrast=c(-.5,.5,-.5,.5))
```
To look at the genes with the lowest FDR-corrected p-values:
```
topTags(lrtTrt)
```
And to write all of the test results to file
```
write.table(topTags(lrtTrt, n=16877), 'edgeR_TrtLRT.txt',
            sep='\t')
```
How many genes are differentially expressed by Treatment?

Now create a smear plot for the Treatment comparison, but make sure the axis limits are the same to make comparison easier.

How is it different from the Population smear plot? What does this say about the relative effects of stickleback population and microbiota treatment on gut gene expression in this experiment?

# 4. Testing effects of Population and Treatment on differential expression using the general linear model approach in limma (with voom()) (20 points total)
Note that limma should already have been installed with edgeR.

Repeat everything you did for the edgeR analysis up to the dispersion estimates, but name the appropriate objects with a “.lv” suffix instead of “.er”

Now use voom to generate “precision weights” from the mean-variance relationship.
```
v.lv <- voom(dge.lv, design.lv)
```
Fit the full general linear model, which includes the voom() precision weights
```
fit.lv <- lmFit(v.lv, design.lv)
fit.lv <- eBayes(fit.lv)
```
Now we perform hypothesis tests for factor level contrasts of interest Population contrast To test for the Population effect we contrast the average of Bt.CV and Bt.GF vs. the average of RS.CV and RS.GF
```
fit.lv.Pop <- contrasts.fit(fit.lv,
                            contrast=c(-.5,-.5,.5,.5))
fit.lv.Pop <- eBayes(fit.lv.Pop)
```
To look at the genes with the lowest FDR-corrected p-values:
```
topTable(fit.lv.Pop)
```
And to write all of the test results to file
```
write.table(topTable(fit.lv.Pop, n=16877),
                     'limma_PopLRT.txt', sep='\t')
```
Summarize how many genes are differentially expressed by Population, using an FDR = 0.05. In this case, practice subsetting topTable(fit.lv.Pop, n=16877) using the appropriate logical operators

Treatment contrast To test for the Treatment effect we contrast the average of Bt.CV and RS.CV vs. the average of Bt.GF and RS.GF
```
fit.lv.Trt <- contrasts.fit(fit.lv,
                            contrast=c(-.5,.5,-.5,.5))
fit.lv.Trt <- eBayes(fit.lv.Trt)
```
To look at the genes with the lowest FDR-corrected p-values:

```
topTable(fit.lv.Trt)
```
And to write all of the test results to file:
```
write.table(topTable(fit.lv.Trt, n=16877),
                     'limma_TrtLRT.txt', sep='\t')
```
Summarize how many genes are differentially expressed by Treatment, using an FDR = 0.05. Again, subset topTable(fit.lv.Trt, n=16877)using the appropriate logical operators

# 5. Compare the overlap in differentially expressed genes between edgeR and limma-voom (10 points total)

Make 2 Venn Diagrams to show the overlap of sig. genes (FDR=.05) between edgeR and limma approaches. One diagram should characterize the Population test, and one the Treatment test. You will need to find the intersection of genes between the two lists and use a package like “VennDiagram”.

Briefly comment on any differences in the numbers of differentially expressed genes between the edgeR and limma-voom approaches. In general do the approaches agree? Is one substanitally more conservative than the other?

