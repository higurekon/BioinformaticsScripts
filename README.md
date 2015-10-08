README

These are all bioinformatics scripts that will work if you run them in the Workspace directory.
parse_tcga_files.sh must be run with a large file 10+ GB file so I have not included that in here, but will try to include a more workable example in the near future.
Input files used to run this script are very large originally, and have been truncated to a smaller example.
These scripts were made to run multiple cancer types, although I am only using colorectal cancer as an example here.

snp_comparator is a script that compares the genomes of two strains of E. Coli for mutations and uses a mutation score per gene based on its SNP count, and maps that onto difference in metabolite concentration.

immune_compendium takes a RNAseq dataset and performs a hierarchical cluster and heatmap on it to illuminate clusters of coexpressing immune-related genes.

combined_score_test and quantile_cutoff take expression levels of RORgT and target genes and then use cutoffs to partition the data into quadrants (high vs high, high vs low, low vs high, low vs low). A kaplain-meier survival curve is drawn to show the prognosis of each quadrant of patient. quantile_cutoff tests each target gene individually, while combined_score_test uses a custom score based on all the target genes.

correlation_of_cancer_with_rorgt looks at the correlation of RORgT and its target genes in a cancer type.

correlation_cutoff uses a cutoff that is based on a regression line between RORgT expression and IL17A expression, the best known target of RORgT. The rationale is that data that is closer to this regression line are patients that more likely have active RORgT. A survival curve is plotted for this analysis.