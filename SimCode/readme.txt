In this folder you'll find three R programs that replicate each of the simulation investigations found in Section 4 of "Flexible Bayesian survival modeling with semiparametric time-dependent and shape-restricted covariate effects."

The programs are named after the Section for with they correspond section, and are further annotated:
SimCode4.1.R: Replicates the simulation in Section 4.1
SimCode4.2.R: Replicates the simulation in Section 4.2
SimCode4.3.R: Replicates the simulation in Section 4.3

Each program begins with JAGS code to implement the relevant statistical models, then a simulation function that conducts a single iteration, then code to run these simulation iterations parallelly, and finally code to compile the results.