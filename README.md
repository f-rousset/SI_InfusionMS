### Scripts and results for ["A new iterative framework for simulation-based population genetic inference with improved coverage properties of confidence intervals"](https://www.biorxiv.org/content/10.1101/2024.09.30.615940v3)
By Rousset, Leblois, Estoup & Marin 

The Rmarkdown file [`Infusion.Rmd`](InfusionMS.Rmd) can be used to generate the **Figures and Tables** 
describing the results of this paper, using information saved in subdirectories.  The Rmarkdown file and its [html output](InfusionMS.html), can be consulted to match a given 
subdirectory to a given Table or Figure from the ms.

Raw real **data** used in the study are provided as
[`diy2inf_simuls/Harmonia/data_HA_INFUSION.txt`](diy2inf_simuls/Harmonia/data_HA_INFUSION.txt) and 
[`diy2inf_simuls/admixtOutOfA/human_snp_all22chr_maf1__10inds_per_sample.snp`](diy2inf_simuls/admixtOutOfA/human_snp_all22chr_maf1__10inds_per_sample.snp).
The same directories contain additional files that may be needed to run the simulations, such as 
`statobsRF.txt` that contains the summary statistics for the data.

Scripts that can be used to reproduce the **simulations** are provided in directories
[`toyTests`](toyTests/) and [`diy2inf_simuls`](diy2inf_simuls/). These directories themselves contain nested 
subdirectories. Each terminal subdirectory corresponds to a simulation scenario.
Files from parent directories should be copied in a terminal subdirectory in order 
to reproduce simulations for the corresponding simulation scenario. In particular, 
the two `generic_workflow.R` files are each a master script file whose execution depends 
on the name of the terminal subdirectory where it is copied.

The terminal subdirectories also contain saved summaries of the inferences 
for all simulation replicates for a given simulation scenario, used by `Infusion.Rmd`.
The [`diy2inf`](diy2inf) subdirectory contain definitions of `R` functions used either by the simulation 
code or by the Rmarkdown script.

To run the `Infusion.Rmd` on another computer or directory, there are two directories 
which need to be specified in the first R chunks of the Rmarkdown file:

* `correctdir` : the base directory where the files from this git repository are copied. `correctdir` can easily be changed but the subdirectory architecture of the files from the git repository should not be modified.
* `MSdir` : the directory where figures and tables will be written if the Rmarkdown script is run.

Further, the `latex2rtf` and the `ggplot2` packages should be installed.

Files related to `sbi` simulations were added in revision, in `SNLE/` subdirectories for three simulations scenarios. They include the `python` scripts used to run the simulations, e.g., [`snle_N_7from17.py`](diy2inf_simuls/admixtOutOfA/N_7from17/SNLE/allfilesforeachjob/snle_N_7from17.py), which can be run in an environment defined by
```
conda create -n sbi_env python=3.12 && conda activate sbi_env
pip install torch
pip install numpy
pip install sbi==0.25.0
```

  
  
