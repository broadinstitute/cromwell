Most users should not need to build Cromwell and can use pre-built Cromwell [releases](Getting).

If for some reason you require a non-release version of Cromwell or are developing new Cromwell
features or fixes, the following are required to build Cromwell from source:

* [Scala 2.12](http://www.scala-lang.org/)
* [SBT 1.x](https://www.scala-sbt.org/)
* [Java 8](http://www.oracle.com/technetwork/java/javase/overview/java8-2100321.html)
* [Git](https://git-scm.com/)

You can also use the [development image](https://github.com/broadinstitute/cromwell/tree/develop/scripts/docker-develop), and build a development container to work inside:

```bash
$ docker build -t cromwell-dev Dockerfile
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

`sbt assembly` will build the runnable Cromwell JAR in `server/target/scala-2.12/` with a name like `cromwell-<VERSION>.jar`.

To build a [Docker](https://www.docker.com/) image, run:

```bash
$ sbt server/docker
```

This will build and tag a Docker image with a name like `broadinstitute/cromwell:<VERSION>-SNAP`.

###IntelliJ gotchas

If you develop Cromwell using IntelliJ IDEA, you may see the following error message in IntelliJ:

```
Error occurred during initialization of VM
Initial heap size set to a larger value than the maximum heap size
```

This is caused by IntelliJ overriding the `Xmx` (max heap) option in our custom `.sbtopts` file to a value that's below 
our `Xms` (initial heap) value.  We cannot tell IntelliJ not to do this, but we can set it to a compatible value.
Go to **Preferences** -> **Build, Execution, Deployment** -> **sbt** -> **JVM** and increase the value of 
`Maximum heap size, MB` above our `Xms` value of 2 GB.  A `Maximum heap size` value of `2100 MB` will be sufficient.

