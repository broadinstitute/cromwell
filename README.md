[![Build Status](https://travis-ci.org/broadinstitute/wdl4s.svg?branch=develop)](https://travis-ci.org/broadinstitute/wdl4s?branch=develop)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/wdl4s/badge.svg?branch=develop)](https://coveralls.io/r/broadinstitute/wdl4s?branch=develop)

# Scala binding API for WDL

This repository provides scala tools to parse a [WDL](https://github.com/broadinstitute/wdl) file and transform it into a scala object hierarchy.

## Installation

wdl4s is hosted on The Broad Institute's [Artifactory Repository]()

```
resolvers ++= Seq(
  "Broad Artifactory Releases" at "https://artifactory.broadinstitute.org/artifactory/libs-release/",
  "Broad Artifactory Snapshots" at "https://artifactory.broadinstitute.org/artifactory/libs-snapshot/"
)
```

Add to `libraryDependencies` with:

```
"org.broadinstitute" %% "wdl4s" % "0.3",
```

Or add a snapshot release in the format `<version>-<git-hash8>-SNAPSHOT`:

```
"org.broadinstitute" %% "wdl4s" % "0.3-e1d8072-SNAPSHOT",
```


To use Squants in your Maven project add the following dependency

```xml
<dependency>
    <groupId>org.broadinstitute</groupId>
    <artifactId>wdl4s_2.11</artifactId>
    <version>0.3</version>
</dependency>
```

## Scaladoc

* [0.3](http://broadinstitute.github.io/wdl4s/0.3)

## Usage

The main entry point into the parser is the `WdlNamespace` object.  A [WDL](https://github.com/broadinstitute/wdl) file is considered a namespace, and other namespaces can be included by using the `import` statement (but only with an `as` clause).

Namespaces have within them 

```scala
import java.io.File
import wdlscala.NamespaceWithWorkflow

object main {
  def main(args: Array[String]) {
    val ns = NamespaceWithWorkflow.load("""
    |task a {
    |  command { ps }
    |}
    |workflow wf {
    | call a
    |}""".stripMargin)

    println(s"Workflow: ${ns.workflow.name}")
    ns.workflow.calls foreach {call =>
      println(s"Call: ${call.name}")
    }

    ns.tasks foreach {task =>
      println(s"Task: ${task.name}")
      println(s"Command: ${task.commandTemplate}")
    }
  }
}
```

To access only the parser, use the `AstTools` library, as follows:

```scala
import java.io.File
import wdlscala.AstTools
import wdlscala.AstTools.EnhancedAstNode

object main {
  def main(args: Array[String]) {
    /* Create syntax tree from contents of file */
    val ast = AstTools.getAst(new File(args(0)))

    /* Second parameter is a descriptor about where the first string came from.
     * Most of the time this would be the URI of where the text was loaded from,
     * but there are no restrictions on what the string can be.
     */
    val ast2 = AstTools.getAst("workflow simple {}", "string")

    /* Print the AST */
    println(ast.toPrettyString)

    /* Traverse the tree to find all Task definitions */
    val taskAsts = AstTools.findAsts(ast, "Task") foreach {ast =>
      println(s"Task name: ${ast.getAttribute("name").sourceString}")
    }
  }
}
```

# Developer Zone

This section is intended for developers of wdl4s.

## Generate WDL Parser

If there are any changes to `src/main/resources/grammar.hgr`, then the parser located at `src/main/java/wdl4s/parser/WdlParser.java` needs to be regenerated.  Any changes to the grammar file should result in a regeneration of the parser and then run the unit tests. Changing the AST could be disruptive if keys are renamed or objects restructured too much. It's best to find these issues as soon as possible.

To regenerate the parser, install the latest version of [Hermes](http://github.com/scottfrazer/hermes), then run the following command within this directory:

```
hermes generate src/main/resources/grammar.hgr \
  --language=java \
  --directory=src/main/java \
  --name=wdl \
  --java-package=wdl4s.parser \
  --java-use-apache-commons \
  --java-imports=org.apache.commons.lang3.StringEscapeUtils \
  --header
```
