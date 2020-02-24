**Cromwell Releases**

Cromwell releases are available at the [GitHub Releases](https://github.com/broadinstitute/cromwell/releases/latest) page. 
You are strongly encouraged to use the latest release of Cromwell whenever possible.

Cromwell is distributed as a conda package on [conda-forge](https://conda-forge.org/).
These instructions need to be followed for [installing the miniconda distribution](https://docs.conda.io/en/latest/miniconda.html) and 
[activating the conda-forge channel](https://conda-forge.org/#about). After this Cromwell can be installed in the 
base environment with `conda install cromwell` or a separate environment for Cromwell can be created with 
`conda create -n cromwell cromwell`. If you are using Cromwell for bioinformatics workflows, you might like to take
a look at [bioconda](http://bioconda.github.io)  as well. 
The conda installation of Cromwell comes with a wrapper that locates the jar for you and allows for running Cromwell or Womtool with a 
`cromwell run`, `womtool validate` or other command. Conda also installs the required Java dependency 
in the environment automatically.

Mac users with Homebrew can also get Cromwell with the command `brew install cromwell`.

This documentation frequently refers to a "Cromwell jar" with a name like `cromwell-<version>.jar`. 
This is the main artifact in Cromwell releases that contains all executable Cromwell code and default configuration.   

[Java 8](http://www.oracle.com/technetwork/java/javase/overview/java8-2100321.html) is required to run Cromwell.

For users running a cromwell server [a docker image](https://hub.docker.com/r/broadinstitute/cromwell) has been made available.
