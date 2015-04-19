Cromwell
========

Workflow engine using [WDL](SPEC.md) as the workflow and task language.

Generate WDL Parser
-------------------

Install the latest version of [Hermes](http://github.com/scottfrazer/hermes), then run the following command within this directory:

```
hermes generate src/main/resources/grammar.hgr --language=java --directory=src/main/java --name=wdl --java-package=cromwell.parser --java-use-apache-commons
```

Architecture
------------

The architecture is split into three layers, from bottom to top:

### cromwell.parser

Contains only the WDL parser to convert WDL source code to an abstract syntax tree

### cromwell.binding

Contains code that takes an abstract syntax tree and returns native Scala object representations of those ASTs.  This layer will also have functions for evaluating expressions when support for that is added

### cromwell.engine

Contains the Akka code and actor system to execute a workflow.
