## Running GWAS using regenie nextflow pipeline


The pipeline can run several phenotypes at once, but all phenotype files need to have the same number of phenotypes (`num_phenotypes_per_file` in `nextflow.config`).

Given a file with many phenotypes, one can use the script `utils/split_pheno.sh` to split it into multiple files.
For example,

```sh
cut -f1,2,3-102 phecodeX_tble.tsv > pheno100.tsv
utils/split_pheno.sh pheno100.tsv 20
```

extracts the first 100 phenotypes and splits them into files with 20 phenotypes each.


All data that is to be accessible by the nextflow pipeline needs to be in the workspace bucket.
We, therefore, need to copy the phenotype files into the workspace bucket.



Set parameters in `nextflow.config` and run the pipeline using 

```sh
~/nextflow run main.nf -profile gcb,spot -with-report report.html -with_trace trace.tsv
```

- If the spot machines are interrupted, the pipeline run can be continued using the `-resume` flag.
- To check for errors, look at the `.nextflow.log` file in the run directory.
- To check recent runs, do `~/nextflow log`
- To check certain run, do `~/nextflow log <RUN> -f workdir,name,status`




