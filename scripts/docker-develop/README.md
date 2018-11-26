# Cromwell Docker Development

This is a container intended to help with development of Cromwell, in the
case that you don't want to install the dependencies on your host. It
includes the required software mentioned in [the developer docs](http://cromwell.readthedocs.io/en/develop/Building/) but with additional instructions to interact with the repository.

As with the build instructions, first start by cloning the Cromwell repository from GitHub:

```bash
$ git clone git@github.com:broadinstitute/cromwell.git
```

And change into this directory:

```bash
$ cd cromwell/scripts/docker-develop
```

This is where this README.md sits with a [Dockerfile](Dockerfile)

## Step 1. Build the Container
From this folder with the Dockerfile, build the container. In the command below
we are calling it `cromwell-dev` and you can choose to change this name if you want.

```bash
$ docker build -t cromwell-dev .
```

and change back to the root of the repository

```bash
cd ../../
```

## Step 2. Shell into Working Environment
If you require a specific version of Cromwell as a starting point, do the appropriate `git checkout` now. 

You are next going to want to bind the cromwell
source code to be somewhere in the container. Actually, we have a `/code` directory
ready for you to make this easy. Run this command from the base of the repository:

```bash
$ docker run -v $PWD/:/code -it cromwell-dev bash
```

## Step 4. The Development Steps

You are going to use the command `sbt assembly` to create a runnable Cromwell JAR. The output
will be `server/target`. Remember, here we are inside the container with scala and sbt:


```bash
$ sbt assembly
```

`sbt assembly` will build the runnable Cromwell JAR in `server/target/scala-2.12/` with a name like `cromwell-<VERSION>.jar`.

You can then interact with it in the container, or on your host if you like. Remember that
you have Java already in the container, so it makes sense to develop there.

