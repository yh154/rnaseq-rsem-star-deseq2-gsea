fun_run_gsea <- function(rnk, jar, gmt, nperm, nplot, smax, smin){
  prefix = gsub(".rnk","",basename(rnk))
  log <- file(sprintf("gsea/%s.log",prefix), open="wt")
  sink(log)
  sink(log, type="message")
  gsea.cmd <- sprintf("java -Xmx8g -cp %s xtools.gsea.GseaPreranked -gmx %s -collapse false -rnk %s -nperm %i -scoring_scheme weighted -rpt_label %s -make_sets true -plot_top_x %i -rnd_seed 521 -set_max %i -set_min %i -zip_report true -out gsea/ -gui false > gsea/%s.log",
  jar, gmt, rnk, nperm, prefix, nplot, smax, smin, prefix)
  system(gsea.cmd)
}

dir.create("gsea")

jar <- snakemake@params[["jar"]]
gmt <- snakemake@params[["gmt"]]
rnk <- snakemake@input

# numbers of permutations (should use 1000 or 10000)
nperm <- snakemake@params[["no_perm"]]
nplot <- snakemake@params[["no_plot"]]
smax <- snakemake@params[["set_max"]]
smin <- snakemake@params[["set_min"]]
seed <- snakemake@params[["seed"]]


## build GSEA command
sapply(rnk, fun_run_gsea, jar=jar, gmt=gmt, nperm=nperm, nplot=nplot, smax=smax, smin=smin )
