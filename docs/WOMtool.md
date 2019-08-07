Command line utilities for interacting with the Workflow Object Model (WOM). You can download the latest WOMtool from the [Cromwell releases page on Github](https://github.com/broadinstitute/cromwell/releases/latest).

## Requirements

The following is the toolchain used for development of womtool.  Other versions may work, but these are recommended.

* [Scala 2.12](http://www.scala-lang.org/)
* [SBT 1.x](https://www.scala-sbt.org/)
* [Java 8](http://www.oracle.com/technetwork/java/javase/overview/java8-2100321.html)
* [Git](https://git-scm.com/)

## Building

`sbt assembly` will build a runnable JAR in `target/scala-2.11/`

Tests are run via `sbt test`.  Note that the tests do require Docker to be running.  To test this out while downloading the Ubuntu image that is required for tests, run `docker pull ubuntu:latest` prior to running `sbt test`

## Command Line Usage

Run the JAR file with no arguments to get the usage message:

`$ java -jar womtool.jar`

```bash
java -jar /path/to/womtool.jar <action> <parameters>

Actions:
validate [--list-dependencies] <WDL file>

  Performs full validation of the WDL file including syntax
  and semantic checking. -l or --list-dependencies is an optional flag to 
  list files referenced in import statements.

inputs <WDL file>

  Print a JSON skeleton file of the inputs needed for this
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

graph <WDL file>

  Reads a WDL file against the grammar and prints out a
  .dot of the DAG if it is valid, and a syntax error
  otherwise.

womgraph <WDL or CWL file> [ancillary files]

  Reads a WDL or CWL file from the first argument and
  converts it to a WOM representation then prints out a graph
  of the WOM produced.
  Any imported files can be supplied as subsequent arguments.
```

### `validate`

Given a WDL file, this runs the full syntax checker over the file and resolves imports in the process.  If any syntax errors are found, they are written out.  Otherwise the program exits.

Error if a `call` references a task that doesn't exist:

`$ java -jar womtool.jar validate 2.wdl`

```
ERROR: Call references a task (BADps) that doesn't exist (line 22, col 8)

  call BADps
       ^
```

Error if namespace and task have the same name:

`$ java -jar womtool.jar validate 5.wdl`

```
ERROR: Task and namespace have the same name:

Task defined here (line 3, col 6):

task ps {
     ^

Import statement defined here (line 1, col 20):

import "ps.wdl" as ps
                   ^
```

##### --list-dependencies or -l flag

For a successful validation, this will output the list of files referenced in import statements in workflows and their subworkflows.

`$ java -jar womtool.jar validate -l myWdl.wdl`

```hocon
Success!
List of Workflow dependencies are:
/path/to/my/import/myImport.wdl
/path/to/another/import/anotherImport.wdl
https://path-to-http-import/httpImport.wdl
```

### `inputs`

Examine a WDL file with one workflow in it, compute all the inputs needed for that workflow and output a JSON template that the user can fill in with values.  The keys in this document should remain unchanged.  The values tell you what type the parameter is expecting.  For example, if the value were `Array[String]`, then it's expecting a JSON array of JSON strings, like this: `["string1", "string2", "string3"]`

`$ java -jar womtool.jar inputs 3step.wdl`


```
{
  "three_step.cgrep.pattern": "String"
}
```

This inputs document is used as input to the `run` subcommand.

### `highlight`

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

### `parse`

Given a WDL file input, this does grammar level syntax checks and writes out the resulting abstract syntax tree.

`$ echo "workflow wf {}" | java -jar womtool.jar parse /dev/stdin`


```
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

`$ java -jar womtool.jar highlight test.wdl html`

```
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

### `graph`
 
The syntax of the graph command is:
```
womtool graph [--all] wdlFile.wdl
```

Given a WDL file input, command generates the data-flow graph through the system in `.dot` format.

For example the fork-join WDL:
```
task mkFile {
  command {
    for i in `seq 1 1000`
    do
      echo $i
    done
  }
  output {
    File numbers = stdout()
  }
  runtime {docker: "ubuntu:latest"}
}

task grep {
  String pattern
  File in_file
  command {
    grep '${pattern}' ${in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
  runtime {docker: "ubuntu:latest"}
}

task wc {
  File in_file
  command {
    cat ${in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
  runtime {docker: "ubuntu:latest"}
}

task join {
  Int grepCount
  Int wcCount
  command {
    expr ${wcCount} / ${grepCount}
  }
  output {
    Int proportion = read_int(stdout())
  }
  runtime {docker: "ubuntu:latest"}
}

workflow forkjoin {
  call mkFile
  call grep { input: in_file = mkFile.numbers }
  call wc { input: in_file=mkFile.numbers }
  call join { input: wcCount = wc.count, grepCount = grep.count }
  output {
    join.proportion
  }
}
```

Produces the DAG:
```
digraph forkjoin {
  "call forkjoin.mkFile" -> "call forkjoin.wc"
  "call forkjoin.mkFile" -> "call forkjoin.grep"
  "call forkjoin.wc" -> "call forkjoin.join"
  "call forkjoin.grep" -> "call forkjoin.join"
}
```

#### The `--all` flag

If this flag is set, all WDL graph nodes become nodes in the generated DAG, even if they are not "executed". Typically this will mean task declarations and call outputs. 
For example in the above example, with `--all` you would get:

```
digraph forkjoin {
  "call forkjoin.grep" -> "String forkjoin.grep.pattern"
  "call forkjoin.grep" -> "output { forkjoin.grep.count = read_int(stdout()) }"
  "call forkjoin.grep" -> "File forkjoin.grep.in_file"
  "call forkjoin.wc" -> "output { forkjoin.wc.count = read_int(stdout()) }"
  "call forkjoin.grep" -> "call forkjoin.join"
  "call forkjoin.wc" -> "File forkjoin.wc.in_file"
  "call forkjoin.mkFile" -> "call forkjoin.grep"
  "call forkjoin.join" -> "output { forkjoin.join.proportion = read_int(stdout()) }"
  "call forkjoin.join" -> "Int forkjoin.join.wcCount"
  "call forkjoin.wc" -> "call forkjoin.join"
  "call forkjoin.mkFile" -> "output { forkjoin.mkFile.numbers = stdout() }"
  "call forkjoin.mkFile" -> "call forkjoin.wc"
  "call forkjoin.join" -> "Int forkjoin.join.grepCount"
}
```
