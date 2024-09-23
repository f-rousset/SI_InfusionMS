Scripts and results for "Better confidence intervals in simulation-based inference"
By Rousset, Leblois, Estoup & Marin

Scripts that can be used to reproduce the simulations are provides in directories
toyTests/ and diy2inf_simuls/ . These directories themselves contain nested 
subdirectories. Each terminal subsdirectory corresponds to a a simulation scenario.
Files from parent directories should be copied in a terminal subdirectory in order 
to reproduce simulations for the corresponding simulation scenario. In particular, 
the two 'generic_workflow.R' files are each a master script file whose execution depends 
on the name of the terminal subdirectory where it is copied.

The terminal subdirectories also contain saved summaries of the inferences 
for all simulation replicates for a given simulation scenario.
The Rmarkdown file Infusion.Rmd can be used to generate the Figures and Tables 
describing the results, using these saved summaries.  
The Rmarkdown file and its html output, should be consulted to match a given 
subdirectory to a given Table or Figure.

The diy2inf/ subdirectory contain definitions of R functions used either by the simulation 
code or by the Rmarkdown script.
  
  
