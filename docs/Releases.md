# Cromwell Releases

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

For users running a Cromwell server [a docker image](https://hub.docker.com/r/broadinstitute/cromwell) has been made available.

### Apple Silicon support statement (updated 2020-11-17)

#### Cromwell JAR works out of the box

The Cromwell JAR works on any standard Java installation. A user can install an x86 Java runtime on an Apple Silicon Mac and the Rosetta 2 translation layer runs Cromwell at near-native speed.

Once natively-compiled Java runtimes become available, performance will increase with no change in functionality. 

#### Docker Desktop support is in progress

The Cromwell Docker image will not run on M1 Macs until Docker Desktop ships the appropriate update. For more details, please see [their official announcement](https://www.docker.com/blog/apple-silicon-m1-chips-and-docker/).

By extension, the absence of Docker means that Cromwell's local Docker backend is not yet supported.

Even when Docker Desktop goes native on Apple Silicon, any tool images running on the local backend will need to cross-compile for the x86 and Arm architectures. This is because the Rosetta 2 translation layer [does not support virtualization](https://developer.apple.com/documentation/apple_silicon/about_the_rosetta_translation_environment). Please contact the tool maintainers for more information. 
