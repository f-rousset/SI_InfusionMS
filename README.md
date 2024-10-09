### Scripts and results for ["Better confidence intervals in simulation-based inference"](https://www.biorxiv.org/content/10.1101/2024.09.30.615940v1)
By Rousset, Leblois, Estoup & Marin 

The Rmarkdown file [`Infusion.Rmd`](InfusionMS.Rmd) can be used to generate the **Figures and Tables** 
describing the results of this paper, using information saved in subdirectories.  The Rmarkdown file and its [html output](InfusionMS.html), can be consulted to match a given 
subdirectory to a given Table or Figure from the ms.

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
* `correctdir` : the base directory where the files from this git repository are copied. `correctdir` can easily be changed but the subdirectory architecture of 
the files from the git repository should not be modified.
* `MSdir` : the directory where figures and tables will be written if the Rmarkdown script is run.

  
  
