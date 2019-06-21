Most users should not need to build Cromwell and can use pre-built Cromwell [releases](Getting).

If for some reason you require a non-release version of Cromwell or are developing new Cromwell
features or fixes, the following are required to build Cromwell from source:

* [Scala 2.12](http://www.scala-lang.org/)
* [SBT 1.x](https://www.scala-sbt.org/)
* [Java 8](http://www.oracle.com/technetwork/java/javase/overview/java8-2100321.html)
* [Git](https://git-scm.com/)

You can also use the [development image](https://github.com/broadinstitute/cromwell/tree/develop/scripts/docker-develop), and build a development container to work inside:

```bash
$ docker build -t cromwell-dev .
$ docker run -it cromwell-dev bash
```

First start by cloning the Cromwell repository from GitHub:

```bash
$ git clone git@github.com:broadinstitute/cromwell.git
```

Next change into the `cromwell` directory:

```bash
$ cd cromwell
```

If you require a specific version of Cromwell as a starting point, do the appropriate `git checkout` now. 

Finally build the Cromwell jar:

```bash
$ sbt assembly
```

NOTE: This command will run for a long time the first time.  
NOTE: Compiling will not succeed on directories encrypted with ecryptfs (ubuntu encrypted home dirs for example), due to long file paths.

`sbt assembly` will build the runnable Cromwell JAR in `server/target/scala-2.12/` with a name like `cromwell-<VERSION>.jar`.

To build a [Docker](https://www.docker.com/) image, run:

```bash
$ sbt server/docker
```

This will build and tag a Docker image with a name like `broadinstitute/cromwell:<VERSION>-SNAP`.

