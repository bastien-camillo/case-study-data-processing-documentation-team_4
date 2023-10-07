First, we activated our environment and installed some R packages using mamba​

```bash
conda activate day3
mamba install r-tidyverse r-ggrepel
```

Then we set up a working folder and copied the data 

```bash
mkdir popgen
cd popgen
cp ~/course/data/day3/popgen/* .
```

Our dataset consisted of three files

```bash
modern_polar_mexican.bed
modern_polar_mexican.bim
modern_polar_mexican.fam
```

We used PLINK for SNP genotype data managment and analysis. 

```bash
plink --bfile modern_polar_mexican --missing --out modern_polar_mexican
```

To run PCA we did this

```bash
smartpca -p modern_all.smartpca.par | tee modern_all.smartpca.log
```

With the output of parameter file:

```bash
genotypename:    modern_polar_mexican.bed
snpname:         modern_polar_mexican.bim
indivname:       modern_polar_mexican.fam
evecoutname:     modern_all.evec
evaloutname:     modern_all.eval
familynames:     NO
numoutevec:      20
numthreads:	 2
numoutlieriter:	 0
poplistname:   	 modern_all.pops.txt
lsqproject:  YES
pordercheck:  NO
autoshrink: YES
```

We then used R to visualize our results

```r
#!/usr/bin/Rscript --vanilla

#########################################
## A simple R script to plot smartpca  ##
#########################################

## --------------------------------------------------------
## libraries

library(tidyverse)
library(ggrepel)


## --------------------------------------------------------
## arguments

args <- commandArgs(trailingOnly = TRUE)
evec_file <- args[1]
label_file <- args[2]
out_file <- args[3]


## --------------------------------------------------------
## read and plot the data

## eigenvector data table and helpers
d <- read_table(
    file = evec_file,
    skip = 1,
    col_names = FALSE
)

samples_label <- scan(
    file = label_file,
    what = character()
)

## some reformatting
n_pcs <- ncol(d) - 2
colnames(d)[1] <- "sample_id"
colnames(d)[ncol(d)] <- "pop_id"
colnames(d)[2:(ncol(d) - 1)] <- paste("PC", 1:n_pcs, sep = "")

d_bg <- d %>%
    filter(!sample_id %in% samples_label)

d_label <- d %>%
    filter(sample_id %in% samples_label)

pcs <- matrix(paste("PC", 1:n_pcs, sep = ""), nrow = 2)

## plot
pdf(
    file = out_file,
    width = 6,
    height = 5
)

walk(seq_len(ncol(pcs)), ~ {
    p <- ggplot(d_bg, aes(
        x = !!sym(pcs[1, .x]),
        y = !!sym(pcs[2, .x])
    ))

    print(p +
        geom_point(
            aes(
                color = pop_id,
                fill = pop_id
            ),
            size = 1.5,
            alpha = 0.9
        ) +
        geom_point(
            size = 2,
            alpha = 1,
            fill = "grey40",
            shape = 21,
            data = d_label
        ) +
        geom_text_repel(
            aes(
                label = sample_id
            ),
            size = 1.5,
            segment.size = 0.25,
            segment.color = "grey",
            data = d_label,
            max.overlaps = Inf
        ) +
        scale_color_brewer(
            name = "Population",
            palette = "Set2"
        ) +
        scale_fill_brewer(
            name = "Population",
            palette = "Set2"
        ) +
        theme_bw())
})
dev.off()
```


Then another R script to run f statistics

```r
#!/usr/bin/Rscript --vanilla

#################################
## f-statistic with admixtools ##
#################################

## --------------------------------------------------------
## libraries

## here we load the libraries required for the analyses

library(tidyverse)
library(admixtools)
library(ape)
library(ggtree)


## --------------------------------------------------------
## helpers

## we set up a few helper objects for our analysis

populations <- c("Kenai", "SEAK", "West", "Southwest", "East")
outgroup <- "Polar"
populations_all <- c(outgroup, populations, "Mexican")

## here we set the path and prefix of the PLINK format file containing the genotype data to analyse
## in order to use admixtools on PLINK files, 
## the population labels have to be specified in the 'family ID' column (column 1) of the dataset fam file

data_prefix <- "modern_polar_mexican" 


## --------------------------------------------------------
## f4 statistics vs east

## In our first analysis, we want to test whether our ancient mexican bear population (`pop2`)
## is significantly more closely related to the "East" modern population (`pop4`) than to the
## remaining modern populations (vector `pop3`), using the f4 statistic.
## We use the polar bear population as outgroup (`pop1`)
## The `blg_size` parameter sets the size of genomic blocks for the jacknife standard error calculation

r <- qpdstat(
    data = data_prefix,
    pop1 = outgroup,
    pop2 = "Mexican",
    pop3 = populations[populations != "East"],
    pop4 = "East",
    blgsize = 500
)

## plot the results (estimate of f4 and +/- 3 standard errors) for each test population 

pdf(
    file = "f4_east.pdf",
    width = 5,
    height = 2
)
p <- ggplot(r, aes(
    x = est,
    y = pop3
))
p +
    geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.25
    ) +
    geom_errorbarh(
        aes(
            xmin = est - 3 * se,
            xmax = est + 3 * se,
        ),
        height = 0.1,
        linewidth = 0.25,
        color = "grey40"
    ) +
    geom_point(aes(color = pop3),
        size = 2,
    ) +
    scale_color_brewer(
        name = "Population",
        palette = "Set2"
    ) +
    xlab(expression(italic(f)[4])) +
    ylab("Population") +
    theme_bw()
dev.off()

## Examine the output plot. 
## Are there any modern populations with an f4 significantly different from f4=0?
## If so, which ones are they, and how can we interpret this signal?


## --------------------------------------------------------
## admixture f3 for SEAK as ancient west vs east

## In the next analysis, we want to test for evidence of admixture 
## between the two deeply diverged ancestral eastern and western bear clades 
## in the modern "SEAK" population (`pop1`).
## We use the "Southwest" population as representative of the ancestral western clade (`pop3`)
## and test three populations related to the ancestral eastern clade as other sources (`pop2`)

r <- qp3pop(
    data = data_prefix,
    pop1 = "SEAK",
    pop2 = c("Mexican", "East", "Kenai"),
    pop3 = "Southwest",
    blgsize = 500
)

## plot the results (estimate of f3 and +/- 3 standard errors) for each eastern source population 

pdf(
    file = "f3admx_seak.pdf",
    width = 4,
    height = 3
)
p <- ggplot(r, aes(
    x = pop2,
    y = est
))
p +
    geom_hline(
        yintercept = 0,
        linetype = "dashed",
        linewidth = 0.25
    ) +
    geom_errorbar(
        aes(
            ymin = est - 3 * se,
            ymax = est + 3 * se,
        ),
        width = 0.1,
        linewidth = 0.25,
        color = "grey40"
    ) +
    geom_point(
        aes(color = pop2),
        size = 2,
    ) +
    scale_color_brewer(
        name = "Population",
        palette = "Set2"
    ) +
    xlab("Population") +
    ylab(expression(italic(f)[3])) +
    theme_bw()
dev.off()

## Examine the output plot. 
## Are there any eastern source populations with an f3 significantly below zero?
## Can we conclude whether there is evidence for admixture in the "SEAK" population?


## --------------------------------------------------------
## pairwise FST for all populations

## In the final analysis, we calculate FST between pairs of populations
## FST is a classic population genetic statistic measuring the genetic differentiation between populations
## We use all populations for both `pop1` and `pop2` to obtain all pairwise values

r <- fst(
    data = data_prefix,
    pop1 = populations_all,
    pop2 = populations_all,
    adjust_pseudohaploid = TRUE,
    blg_size = 500
)

## here we re-format the results tibble into a matrix of pairwise FST
fst_matrix <- r %>%
    select(-se) %>%
    pivot_wider(
        names_from = "pop2",
        values_from = "est"
    ) %>%
    column_to_rownames("pop1") %>%
    as.matrix()

## finally, we build a neighbor joining tree of all populations based on the FST distances
## and root it using the polar bear population as outgroup
tr <- nj(fst_matrix) %>%
    root(outgroup = "Polar")

## we use the `ggtree` package to visaulize the resulting tree
pdf(
    file = "fst_tree.pdf",
    width = 4,
    height = 3
)
ggtree(tr) +
    geom_tippoint(aes(color = label),
        size = 3
    ) +
    geom_tiplab(
        size = 3,
        offset = 0.02
    ) +
    scale_color_brewer(
        name = "Population",
        palette = "Set2"
    ) +
    xlim(c(0, 0.8)) +
    theme_tree()
dev.off()

## Examine the output plot. 
## Which population is closest to the ancient Mexican bear population?
## Compare this tree to the admixture graph in Fig. 1D in Pedersen et al 2021
## How do they differ, and what could be the interpretation?
```
