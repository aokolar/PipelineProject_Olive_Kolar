#load package
library(sleuth)

#read in table created with details from kallisto output
stab = read.table("sleuth_table.txt",header=TRUE)

#initialize sleuth object with sleuth_prep
so = sleuth_prep(stab)

#fit model to compare two conditions (2dpi vs 6dpi)
so = sleuth_fit(so, ~condition, 'full')

#fit reduced model to compare in likelihood test
so = sleuth_fit(so, ~1, 'reduced')

#likelihood ratio test for differential expression between 2dpi / 6dpi
so = sleuth_lrt(so, 'reduced', 'full')

#load dplyr package
library(dplyr)

#pull out test results from sleuth
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

#filter significant results (FDR < .05) and sort by pval
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05) |> dplyr::arrange(pval) 

#write significant results to outfile "fdr05_results"
write.table(sleuth_significant, file="fdr05_results.txt",quote = FALSE,row.names = FALSE)