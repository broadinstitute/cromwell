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

Add the following to `libraryDependencies`:

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

* [0.4](http://broadinstitute.github.io/wdl4s/0.4)

## Usage

All examples are located in `src/main/scala/wdl4s/examples` and particular examples can be run via `sbt`:

```
$ sbt "run-main wdl4s.examples.ex1"
```

### Loading WDL Code

The main entry point into the parser is the `WdlNamespace` object.  A [WDL](https://github.com/broadinstitute/wdl) file is considered a namespace, and other namespaces can be included by using the `import` statement (but only with an `as` clause).

the [WdlNamespace](http://broadinstitute.github.io/wdl4s/0.3/#wdl4s.WdlNamespace$) object has a few `load()` functions for turning WDL source into `WdlNamespace` objects.

If the workflow being loaded contains a `workflow` definition, then the `load()` function will return a [NamespaceWithWorkflow](http://broadinstitute.github.io/wdl4s/0.3/#wdl4s.NamespaceWithWorkflow) and otherwise it will return a [NamespaceWithoutWorkflow](http://broadinstitute.github.io/wdl4s/0.3/#wdl4s.NamespaceWithoutWorkflow).

Example `src/main/scala/wdl4s/examples/ex1.scala`

```scala
val wdl = """
  |task a {
  |  command { ps }
  |}
  |workflow wf {
  | call a
  |}""".stripMargin

val ns = NamespaceWithWorkflow.load(wdl)

println(s"Workflow: ${ns.workflow.unqualifiedName}")
ns.workflow.calls foreach {call =>
  println(s"Call: ${call.unqualifiedName}")
}

ns.tasks foreach {task =>
  println(s"Task: ${task.name}")
  println(s"Command: ${task.commandTemplate}")
}
```

## Using An Import Resolver

WDL code can have `import` statements but wdl4s does not know how to resolve these import statements into WDL source code.

When using `WdlNamespace.load()`, one can pass an optional import resolver which is a `String => WdlSource` function (`WdlSource` is a type alias for `String`).  If the import resolver cannot resolve the import string to WDL source, then it is expected to throw an exception.

Example `src/main/scala/wdl4s/examples/ex2.scala`

```scala
val wdl = """
  |import "some_string"
  |task a {
  |  command { ps }
  |}
  |workflow wf {
  | call a
  |}""".stripMargin

def resolver(importString: String): WdlSource = {
  importString match {
	case "some_string" => "task imported { command {ps} }"
	case s if s.startsWith("http://") =>
	  // issue HTTP request
	  throw new NotImplementedError("not implemented")
  }
}

val ns = NamespaceWithWorkflow.load(wdl, resolver)

ns.tasks foreach {task =>
  println(s"Task: ${task.name}")
}
```

Since the resolver is set up to resolve `some_string` to some static WDL code (`task imported { command {ps} }`, the output of this program will show two tasks in this namespace:

```
Task: a
Task: imported
```

WDL also supports `import "something" as namespace_name` format for imports.  In this case, a sub-namespace will be created where the tasks live

Example `src/main/scala/wdl4s/examples/ex3.scala`

```scala
val wdl = """
  |import "some_string" as my_namespace
  |task a {
  |  command { ps }
  |}
  |workflow wf {
  | call a
  |}""".stripMargin

def resolver(importString: String): WdlSource = {
  importString match {
	case "some_string" => "task imported { command {ps} }"
	case _ => throw new NotImplementedError()
  }
}

val ns = NamespaceWithWorkflow.load(wdl, resolver)

ns.tasks foreach {task =>
  println(s"Task: ${task.name}")
}

ns.namespaces foreach { n =>
  n.tasks.foreach { t =>
	println(s"Imported Task: ${t.name} (from ${n.importedAs.get})")
  }
}
```

Since the WDL `import` statement now has an `as` clause, the top-level namespace only has one task.  The top level namespace also has a sub-namespace called `my_namespace` which has one task.  The output of the program will be:

```
Task: a
Imported Task: imported (from my_namespace)
```

## Resolving Fully-Qualified Names

`WdlNamespace` has a `resolve` method which takes a fully-qualified name string and returns back the object that it refers to.

Example `src/main/scala/wdl4s/examples/ex4.scala`

```scala
val wdl = """
  |task a {
  |  command { ps }
  |}
  |workflow wf {
  | call a
  | call a as b
  |}""".stripMargin

val ns = NamespaceWithWorkflow.load(wdl)

println(ns.resolve("wf.a")) // resolves to Call object for `call a`
println(ns.resolve("wf.b")) // resolves to Call object for `call a as b`
println(ns.findTask("a")) // resolves to Task object for `task a`
```

## Getting Dependencies

`Call` objects can have prerequisites: other `Call`s that need to be completed in order to start executing the `Call` in question.

Example `src/main/scala/wdl4s/examples/ex5.scala`

```scala
val wdl = """
  |task a {
  |  command { ps }
  |  output { File procs = stdout() }
  |}
  |
  |task b {
  |  File s
  |  command { wc -l ${s} }
  |}
  |
  |workflow wf {
  | call a
  | call b {input: s=a.procs}
  |}""".stripMargin

val ns = NamespaceWithWorkflow.load(wdl)

Seq(ns.resolve("wf.a"), ns.resolve("wf.b")) foreach { scope =>
  scope match {
	case Some(c: Call) => println(s"Call '${c.fullyQualifiedName}' prerequisites: ${c.prerequisiteScopes}")
  }
}
```

Since `call b` depends on `call a`, the set of prerequisites for `b` has one element and `a` has zero elements as seen from the standard output:

```
Call 'wf.a' prerequisites: Set()
Call 'wf.b' prerequisites: Set([Call name=a, task=[Task name=a commandTemplate=Vector( ps )}]])
```

## Evaluating Expressions

WDL has its own expression language.  The unevaluated expressions are stored in a `WdlExpression` object, which has an `evaluate()` method.

The `evaluate()` method takes two parameters: a lookup function (`String => WdlValue`) and an implementation of the standard library functions (see `trait WdlStandardLibraryFunctions`)

The lookup function is called for each variable that is encountered during expression evaluation.  The corresponding method in the standard library functions implementation is called for function invocations

Example `src/main/scala/wdl4s/examples/ex6.scala`

```scala
val wdl = """
  |workflow wf {
  |  String a = "foo" + "bar"
  |  String b = "hello " + variable
  |  String c = "hello " + other_variable
  |}""".stripMargin

val ns = NamespaceWithWorkflow.load(wdl)
def lookup(name: String): WdlValue = {
  name match {
	case "variable" => WdlString("world")
	case _ => throw new NoSuchElementException
  }
}
ns.workflow.declarations foreach { decl =>
  val value = decl.expression.get.evaluate(lookup, NoFunctions)
  println(s"Declaration '${decl.toWdlString}' evaluates to: ${value}")
}
```

This evaluates each of the declarations in the workflow.  The last one fails because `other_variable` is not resolved with the `lookup` function that we defined

```
Declaration 'String a = "foo" + "bar"' evaluates to: Success(WdlString(foobar))
Declaration 'String b = "hello " + variable' evaluates to: Success(WdlString(hello world))
Declaration 'String c = "hello " + other_variable' evaluates to: Failure(java.util.NoSuchElementException)
```

## Instantiating Commands

Each `Task`'s command needs to be instantiated from the abstract form in the WDL file to a concrete form.  That means that each expression inside of `${...}` blocks needs to be evaluated.

Example `src/main/scala/wdl4s/examples/ex7.scala`

```scala
val wdl = """
  |task a {
  |  String prefix
  |  Array[Int] ints
  |  command {
  |    python script.py ${write_lines(ints)} > ${prefix + ".out"}
  |  }
  |}
  |workflow wf {
  |  call a
  |}""".stripMargin

val ns = NamespaceWithWorkflow.load(wdl)
val inputs = Map(
  "prefix" -> WdlString("some_prefix"),
  "ints" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(1,2,3,4,5).map(WdlInteger(_)))
)

class CustomFunctions extends WdlFunctions[WdlValue] {
  def write_lines(params: Seq[Try[WdlValue]]): Try[WdlValue] = {
	// Validate `params`, write the result to a file, return file path
	Success(WdlFile("/tmp/array.txt"))
  }
}

ns.findTask("a") foreach { task =>
  println(task.instantiateCommand(inputs, new CustomFunctions).get)
}
```

This will produce the following output:

```
python script.py /tmp/array.txt > some_prefix.out
```

## Accessing the WDL Parser Directly

To access only the parser, use the `AstTools` library, as follows:

Example `src/main/scala/wdl4s/examples/ex8.scala`

```scala
import java.io.File
import wdlscala.AstTools
import wdlscala.AstTools.EnhancedAstNode

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
