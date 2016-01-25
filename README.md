Command line utilities for interacting with WDL

<!---toc start-->

* [Requirements](#requirements)
* [Building](#building)
* [Command Line Usage](#command-line-usage)
  * [validate](#validate)
  * [inputs](#inputs)
  * [highlight](#highlight)
  * [parse](#parse)
* [Getting Started with WDL](#getting-started-with-wdl)

<!---toc end-->

# Requirements

The following is the toolchain used for development of wdltool.  Other versions may work, but these are recommended.

* [Scala 2.11.7](http://www.scala-lang.org/news/2.11.7)
* [SBT 0.13.8](https://github.com/sbt/sbt/releases/tag/v0.13.8)
* [Java 8](http://www.oracle.com/technetwork/java/javase/overview/java8-2100321.html)

# Building

`sbt assembly` will build a runnable JAR in `target/scala-2.11/`

Tests are run via `sbt test`.  Note that the tests do require Docker to be running.  To test this out while downloading the Ubuntu image that is required for tests, run `docker pull ubuntu:latest` prior to running `sbt test`

# Command Line Usage

Run the JAR file with no arguments to get the usage message:

```
$ java -jar wdltool.jar
java -jar wdltool.jar <action> <parameters>

Actions:
validate <WDL file>

  Performs full validation of the WDL file including syntax
  and semantic checking

inputs <WDL file>

  Write a JSON skeleton file of the inputs needed for this
  workflow.  Fill in the values in this JSON document and
  pass it in to the 'run' subcommand.

highlight <WDL file> <html|console>

  Reformats and colorizes/tags a WDL file. The second
  parameter is the output type.  "html" will output the WDL
  file with <span> tags around elements.  "console" mode
  will output colorized text to the terminal
  
parse <WDL file>

Compares a WDL file against the grammar and writes out an
abstract syntax tree if it is valid, and a syntax error
otherwise.  Note that higher-level AST checks are not done
via this sub-command and the 'validate' subcommand should
be used for full validation
```

## validate

Given a WDL file, this runs the full syntax checker over the file and resolves imports in the process.  If any syntax errors are found, they are written out.  Otherwise the program exits.

Error if a `call` references a task that doesn't exist:

```
$ java -jar wdltool.jar validate 2.wdl
ERROR: Call references a task (BADps) that doesn't exist (line 22, col 8)

  call BADps
       ^

```

Error if namespace and task have the same name:

```
$ java -jar wdltool.jar validate 5.wdl
ERROR: Task and namespace have the same name:

Task defined here (line 3, col 6):

task ps {
     ^

Import statement defined here (line 1, col 20):

import "ps.wdl" as ps
                   ^
```

## inputs

Examine a WDL file with one workflow in it, compute all the inputs needed for that workflow and output a JSON template that the user can fill in with values.  The keys in this document should remain unchanged.  The values tell you what type the parameter is expecting.  For example, if the value were `Array[String]`, then it's expecting a JSON array of JSON strings, like this: `["string1", "string2", "string3"]`

```
$ java -jar wdltool.jar inputs 3step.wdl
{
  "three_step.cgrep.pattern": "String"
}
```

This inputs document is used as input to the `run` subcommand.

## highlight

Formats a WDL file and semantically tags it.  This takes a second parameter (`html` or `console`) which determines what the output format will be.

test.wdl
```
task abc {
  String in
  command {
    echo ${in}
  }
  output {
    String out = read_string(stdout())
  }
}

workflow wf {
  call abc
}
```

## parse

Given a WDL file input, this does grammar level syntax checks and writes out the resulting abstract syntax tree.

```
$ echo "workflow wf {}" | java -jar wdltool.jar parse /dev/stdin
(Document:
  imports=[],
  definitions=[
    (Workflow:
      name=<stdin:1:10 identifier "d2Y=">,
      body=[]
    )
  ]
)
```

This WDL file can be formatted in HTML as follows:

```
$ java -jar wdltool.jar highlight test.wdl html
<span class="keyword">task</span> <span class="name">abc</span> {
  <span class="type">String</span> <span class="variable">in</span>
  <span class="section">command</span> {
    <span class="command">echo ${in}</span>
  }
  <span class="section">output</span> {
    <span class="type">String</span> <span class="variable">out</span> = <span class="function">read_string</span>(<span class="function">stdout</span>())
  }
}

<span class="keyword">workflow</span> <span class="name">wf</span> {
  <span class="keyword">call</span> <span class="name">abc</span>
}
```

# Getting Started with WDL

For many examples on how to use WDL see [the WDL site](https://github.com/broadinstitute/wdl/tree/develop#getting-started-with-wdl)
