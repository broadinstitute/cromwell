[![Build Status](https://travis-ci.org/broadinstitute/wdl4s.svg?branch=develop)](https://travis-ci.org/broadinstitute/wdl4s?branch=develop)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/wdl4s/badge.svg?branch=develop)](https://coveralls.io/r/broadinstitute/wdl4s?branch=develop)

# Scala binding API for WDL

The main entry point into the parser is the `WdlNamespace` object.  A [WDL](https://github.com/broadinstitute/wdl) file is considered a namespace, and other namespaces can be included by using the `import` statement (but only with an `as` clause).

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
