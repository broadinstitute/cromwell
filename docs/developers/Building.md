Most users should not need to build Cromwell and can use pre-built Cromwell [releases](Getting).

If for some reason you require a non-release version of Cromwell or are developing new Cromwell
features or fixes, the following are required to build Cromwell from source:

* [Scala 2.13](http://www.scala-lang.org/)
* [SBT 1.x](https://www.scala-sbt.org/)
* [AdoptOpenJDK 11 HotSpot](https://adoptopenjdk.net/)
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

`sbt assembly` will build the runnable Cromwell JAR in `server/target/scala-2.13/` with a name like `cromwell-<VERSION>.jar`. It will also build a runnable Womtool JAR in `womtool/target/scala-2.13/` with a name like `womtool-<VERSION>.jar`.

## Docker

The following Docker build configurations are supported. Most users will want Snapshot, resulting in an image like `broadinstitute/cromwell:<VERSION>-SNAP`.

| Command                                        | Build Type | Debug Tools | Description                          |
|------------------------------------------------|------------|-------------|--------------------------------------|
| `sbt server/docker`                            | Snapshot   | No          | Most common local build              |
| `sbt -Dproject.isDebug=true server/docker`     | Debug      | Yes         | Local build with debugging/profiling |
| `sbt -Dproject.isSnapshot=false server/docker` | Standard   | No          | Reserved for CI: commit on `develop` |
| `sbt -Dproject.isRelease=true server/docker`   | Release    | No          | Reserved for CI: numbered release    |
