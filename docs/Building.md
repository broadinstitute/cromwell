Most users should not need to build Cromwell and can use pre-built Cromwell [releases](Getting).

If for some reason you require a non-release version of Cromwell or are developing new Cromwell
features or fixes, the following are required to build Cromwell from source:

* [Scala 2.12.2](http://www.scala-lang.org/news/2.12.1#scala-212-notes)
* [SBT 1.1.5](https://github.com/sbt/sbt/releases/tag/v1.1.5)
* [Java 8](http://www.oracle.com/technetwork/java/javase/overview/java8-2100321.html)
* [Git](https://git-scm.com/)

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
