Cromwell
========

Workflow engine using [WDL](SPEC.md) as the workflow and task language.

API Documenation
----------------

API Documentation can be found [here](http://broadinstitute.github.io/cromwell/scaladoc).

Generate WDL Parser
-------------------

Install the latest version of [Hermes](http://github.com/scottfrazer/hermes), then run the following command within this directory:

```
hermes generate src/main/resources/grammar.hgr --language=java --directory=src/main/java --name=wdl --java-package=cromwell.parser --java-use-apache-commons
```

Generating and Hosting ScalaDoc
-------------------------------

Essentially run `sbt doc` then commit the generated code into the `gh-pages` branch on this repository

```
$ sbt doc
$ git co gh-pages
$ mv target/scala-2.11/api scaladoc
$ git add scaladoc
$ git commit -m "API Docs"
$ git push origin gh-pages
```

Architecture
------------

![Cromwell Architecture](http://i.imgur.com/kPPTe0l.png)

The architecture is split into four layers, from bottom to top:

### cromwell.parser

Contains only the WDL parser to convert WDL source code to an abstract syntax tree.  Clients should never need to interact with WDL at this level, though nothing specifically precludes that.

### cromwell.binding

Contains code that takes an abstract syntax tree and returns native Scala object representations of those ASTs.  This layer will also have functions for evaluating expressions when support for that is added.

### cromwell.backend

Contains implementations of an interface to launch jobs.  `cromwell.engine` will use this to execute and monitor jobs.

### cromwell.engine

Contains the Akka code and actor system to execute a workflow.  This layer should operate entirely on objects returned from the `cromwell.binding` layer.
