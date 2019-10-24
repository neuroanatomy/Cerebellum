# Cerebellum

![psm_v26_d761_inferior_surface_of_the_cerebellum](https://cloud.githubusercontent.com/assets/2310732/22600049/3762a1ae-ea39-11e6-990f-7cb88209e664.jpg)

## Analysis of cerebellar volumes and areas differences between individuals with ASD and controls

This repository contains the code used for the analyses described in our manuscript "Traut et al. (2017). Cerebellar Volume in Autism: Literature Meta-analysis and Analysis of the Autism Brain Imaging Data Exchange Cohort. Biological Psychiatry. https://doi.org/10.1016/j.biopsych.2017.09.029" (Preprint: https://www.biorxiv.org/content/10.1101/104984v3).

### Meta-analysis of the literature

To compute the meta-analysis and generate the figures, source the file `meta-analysis.R` in a R session.<br>
The studies have to be listed in the file `means-[roi].txt` containing the required fields with the variable `suffix` set to `[roi]`.

### Analysis of Abide

Linear models were fit to evaluate the effect of group and covariables on cerebellum volume using the ABIDE 1 dataset (http://fcon_1000.projects.nitrc.org/indi/abide/).
The script `analyse-abide.R` can be used to run them.

### Meta-analysis of Abide

The script `prepare-meta-abide.R` can be run to prepare Abide data for meta-analysis. It splits data by site and remove individuals to ensure age and sex matching between groups. Mean and standard deviations of volumes for each group by site are then computed and writen to the file `means-abide-Cbl.txt`.
Meta-analysis can then be computed with  `meta-analysis.R` by setting suffix to `abide-Cbl`.

### Other scripts

- `meta-combin.R`: performs a meta-analysis with the literature and abide combined
- `meta-exec.R`: performs a meta-analysis (similar to `meta-analysis.R` but no plot is produced).
- `p-curve.R`: draws p-curves from test statistics. Code was adapted from the source of the p-curve app: `http://p-curve.com/app4/app 4.05.r` by Uri Simonsohn, Leif D. Nelson and Joseph P. Simmons.
- `result-tables-abide.R`: generates report tables for the meta-analysis of the literature.
- `result-tables-mod.R`: generates report tables for the meta-analysis of Abide.
