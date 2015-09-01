[![Build Status](https://travis-ci.org/broadinstitute/cromwell.svg?branch=develop)](https://travis-ci.org/broadinstitute/cromwell?branch=develop)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/cromwell/badge.svg?branch=develop)](https://coveralls.io/r/broadinstitute/cromwell?branch=develop)

Cromwell
========

Workflow engine using [WDL](https://github.com/broadinstitute/wdl/blob/wdl2/SPEC.md) as the workflow and task language.

<!---toc start-->

* [Requirements](#requirements)
* [Building](#building)
* [API Documentation](#api-documentation)
* [Scala API Usage](#scala-api-usage)
* [Command Line Usage](#command-line-usage)
  * [validate](#validate)
  * [inputs](#inputs)
  * [run](#run)
  * [parse](#parse)
  * [server](#server)
* [Getting Started with WDL](#getting-started-with-wdl)
  * [Hello World WDL](#hello-world-wdl)
  * [Modifying Task Outputs](#modifying-task-outputs)
  * [Referencing Files on Disk](#referencing-files-on-disk)
  * [Using String Interpolation](#using-string-interpolation)
  * [Aliasing Calls](#aliasing-calls)
  * [Specifying Inputs and Using Declarations](#specifying-inputs-and-using-declarations)
  * [Using Files as Inputs](#using-files-as-inputs)
  * [Scatter/Gather](#scattergather)
* [Backends](#backends)
  * [Local Filesystem Assumptions / Layout](#local-filesystem-assumptions--layout)
  * [Local Backend](#local-backend)
  * [Sun GridEngine](#sun-gridengine)
  * [Google JES](#google-jes)
* [REST API](#rest-api)
  * [REST API Versions](#rest-api-versions)
  * [POST /workflows/:version](#post-workflowsversion)
  * [GET /workflows/:version/:id/status](#get-workflowsversionidstatus)
  * [GET /workflows/:version/:id/outputs](#get-workflowsversionidoutputs)
  * [GET /workflows/:version/:id/outputs/:call](#get-workflowsversionidoutputscall)
  * [GET /workflows/:version/:id/logs/:call](#get-workflowsversionidlogscall)
  * [GET /workflows/:version/:id/logs](#get-workflowsversionidlogs)
  * [POST /workflows/:version/:id/abort](#post-workflowsversionidabort)
* [Developer](#developer)
  * [Generate WDL Parser](#generate-wdl-parser)
  * [Generating and Hosting ScalaDoc](#generating-and-hosting-scaladoc)
  * [Architecture](#architecture)
    * [cromwell.parser](#cromwellparser)
    * [cromwell.binding](#cromwellbinding)
    * [cromwell.backend](#cromwellbackend)
    * [cromwell.engine](#cromwellengine)

<!---toc end-->

# Requirements

The following is the toolchain used for development of Cromwell.  Other versions may work, but these are recommended.

* [Scala 2.11.7](http://www.scala-lang.org/news/2.11.7)
* [SBT 0.13.8](https://github.com/sbt/sbt/releases/tag/v0.13.8)
* [Java 8](http://www.oracle.com/technetwork/java/javase/overview/java8-2100321.html)

# Building

`sbt assembly` will build a runnable JAR in `target/scala-2.11/`

Tests are run via `sbt test`.  Note that the tests do require Docker to be running.  To test this out while downloading the Ubuntu image that is required for tests, run `docker pull ubuntu:latest` prior to running `sbt test`

# API Documentation

API Documentation can be found [here](http://broadinstitute.github.io/cromwell/scaladoc).

Note that this may not be completely up to date or even useful at this time.

# Scala API Usage

The main entry point into the parser is the `WdlNamespace` object.  A WDL file is considered a namespace, and other namespaces can be included by using the `import` statement (but only with an `as` clause).

```scala
import java.io.File
import cromwell.binding.WdlNamespace

object main {
  def main(args: Array[String]) {
    val ns = WdlNamespace.load(new File(args(0)))
    val ns2 = WdlNamespace.load("workflow wf {}")
    ns.workflows foreach {wf =>
      println(s"Workflow: ${wf.name}")
      wf.calls foreach {call =>
        println(s"Call: ${call.name}")
      }
    }
    ns.tasks foreach {task =>
      println(s"Task: ${task.name}")
      println(s"Command: ${task.command}")
    }
  }
}
```

To access only the parser, use the `AstTools` library, as follows:

```scala
import java.io.File
import cromwell.parser.AstTools
import cromwell.parser.AstTools.EnhancedAstNode

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

# Command Line Usage

Run the JAR file with no arguments to get the usage message:

```
$ java -jar cromwell.jar
java -jar cromwell.jar <action> <parameters>

Actions:

validate <WDL file>

  Performs full validation of the WDL file including syntax
  and semantic checking

inputs <WDL file>

  Print a JSON skeleton file of the inputs needed for this
  workflow.  Fill in the values in this JSON document and
  pass it in to the 'run' subcommand.

run <WDL file> <JSON inputs file>

  Given a WDL file and JSON file containing the value of the
  workflow inputs, this will run the workflow locally and
  print out the outputs in JSON format.

parse <WDL file>

  Compares a WDL file against the grammar and prints out an
  abstract syntax tree if it is valid, and a syntax error
  otherwise.  Note that higher-level AST checks are not done
  via this sub-command and the 'validate' subcommand should
  be used for full validation

highlight <WDL file> <html|console>

  Reformats and colorizes/tags a WDL file. The second
  parameter is the output type.  "html" will output the WDL
  file with <span> tags around elements.  "console" mode
  will output colorized text to the terminal

server

  Starts a web server on port 8000.  See the web server
  documentation for more details about the API endpoints.
```

## validate

Given a WDL file, this runs the full syntax checker over the file and resolves imports in the process.  If any syntax errors are found, they are printed out.  Otherwise the program exits.

Error if a `call` references a task that doesn't exist:

```
$ java -jar cromwell.jar validate 2.wdl
ERROR: Call references a task (BADps) that doesn't exist (line 22, col 8)

  call BADps
       ^

```

Error if namespace and task have the same name:

```
$ java -jar cromwell.jar validate 5.wdl
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
$ java -jar cromwell.jar inputs 3step.wdl
{
  "three_step.cgrep.pattern": "String"
}
```

This inputs document is used as input to the `run` subcommand.

## run

Given a WDL file and a JSON inputs file (see `inputs` subcommand), Run the workflow and print the outputs:

```
$ java -jar cromwell.jar run 3step.wdl inputs.json
[INFO] [06/10/2015 09:20:19.945] [ForkJoinPool-2-worker-13] [akka://cromwell-system/user/$b] Workflow ID: bdcc70e6-e6d7-4483-949b-7c3c2199e26c
[INFO] [06/10/2015 09:20:19.968] [cromwell-system-akka.actor.default-dispatcher-5] [akka://cromwell-system/user/$a/$a] Starting calls: ps
[INFO] [06/10/2015 09:20:20.094] [cromwell-system-akka.actor.default-dispatcher-5] [akka://cromwell-system/user/$a/$a] Starting calls: cgrep, wc
[INFO] [06/10/2015 09:20:20.123] [cromwell-system-akka.actor.default-dispatcher-2] [akka://cromwell-system/user/$b] Workflow complete: Succeeded
{
  "three_step.ps.procs": "/var/folders/kg/c7vgxnn902lc3qvc2z2g81s89xhzdz/T/stdout1272284837004786003.tmp",
  "three_step.cgrep.count": 0,
  "three_step.wc.count": 13
}
```

## parse

Given a WDL file input, this does grammar level syntax checks and prints out the resulting abstract syntax tree.

```
$ echo "workflow wf {}" | java -jar cromwell.jar parse /dev/stdin
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

## highlight

Formats a WDL file and semantically tags it.  This takes a second parameter (`html` or `console`) which determines what the output format will be.

test.wdl
```
task abc {
  String in
  command { echo ${in} }
  output {String out = read_string(stdout())}
}

workflow wf { call abc }
```

This WDL file can be formatted in HTML as follows:

```
$ java -jar cromwell.jar highlight test.wdl html
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

## server

Start a server on port 8000, the API for the server is described in the [REST API](#rest-api) section.

# Getting Started with WDL

If you don't already have a reference to the Cromwell JAR file, compile it with `sbt assembly`, which should produce `cromwell.jar`.

## Hello World WDL

Create a WDL simple file and save it as `hello.wdl`, for example:

```
task hello {
  command {
    echo 'hello ${name}!'
  }
  output {
    File response = stdout()
  }
}

workflow test {
  call hello
}
```

Generate a template `hello.json` file with the `inputs` subcommand:

```
$ java -jar cromwell.jar inputs hello.wdl
{
  "test.hello.name": "String"
}
```

WDL has a concept of fully-qualified names.  In the above output, `test.hello.name` is a fully-qualified name which should be read as: the `name` input on the `hello` call within workflow `test`.  Fully-qualified names are used to unambiguously refer to specific elements of a workflow.  All inputs are specified by fully-qualified names and all outputs are returned as fully-qualified names.

Modify this and save it to `hello.json`:

```
{
  "test.hello.name": "world"
}
```

Then, use the `run` subcommand to run the workflow:

```
$ java -jar cromwell.jar run hello.wdl hello.json
... truncated ...
{
  "test.hello.response": "/Users/sfrazer/projects/cromwell/cromwell-executions/test/c1d15098-bb57-4a0e-bc52-3a8887f7b439/call-hello/stdout8818073565713629828.tmp"
}
```

Since the `hello` task returns a `File`, the result is a file that contains the string "hello world!" in it.

## Modifying Task Outputs

Currently the `hello` task returns a `File` with the greeting in it, but what if we wanted to return a `String` instead?  This can be done by utilizing the `read_string()` function:

```
task hello {
  command {
    echo 'hello ${name}!'
  }
  output {
    String response = read_string(stdout())
  }
}

workflow test {
  call hello
}
```

Now when this is run, we get the string output for `test.hello.response`:

```
$ java -jar cromwell.jar run hello.wdl hello.json
... truncated ...
{
  "test.hello.response": "Hello world!"
}
```

`read_string` is a function in the [standard library](https://github.com/broadinstitute/wdl/blob/wdl2/SPEC.md#standard-library), which provides other useful functions for converting outputs to WDL data types.

## Referencing Files on Disk

So far we've only been dealing with the standard output of a command, but what if it writes a file to disk?  Consider this example:

```
task hello {
  command {
    echo 'hello ${name}!' > out
  }
  output {
    String response = read_string("out")
  }
}

workflow test {
  call hello
}
```

In this example, standard output was redirected to the file "out".  In WDL, `String` data types can be coerced into `File` data types.  Since the `read_string` function needs a file path to read from, it converts `"out"` to a file path and then performs the `read_string` operation.  The output could also be defined as follows:

```
File response = "out"
```

## Using String Interpolation

Sometimes, an output file is named as a function of one of its inputs.

```
task hello {
  command {
    echo 'hello ${name}!' > ${name}.txt
  }
  output {
    String response = read_string("${name}.txt")
  }
}

workflow test {
  call hello
}
```

Here the inputs and outputs are exactly the same as previous examples, however the intermediate output file name of this task is named differently for every invocation.

## Aliasing Calls

Say we wanted to call the `hello` task twice.  Simply adding two `call hello` statements to the body of `workflow test` would result in non-unique fully-qualified names.  To resolve this issue, `call` statements can be aliased using an `as` clause:

```
task hello {
  command {
    echo 'hello ${name}!'
  }
  output {
    String response = read_string(stdout())
  }
}

workflow test {
  call hello
  call hello as hello2
}
```

Now, we need to specify a value for `test.hello2.name` in the hello.json file"

```
{
  "test.hello.name": "world",
  "test.hello2.name": "boston"
}
```

Running this workflow now produces two outputs:

```
$ java -jar cromwell.jar run hello.wdl hello.json
... truncated ...
{
  "test.hello.response": "Hello world!",
  "test.hello2.response": "Hello boston!"
}
```

## Specifying Inputs and Using Declarations

A `call` can have an optional section to define inputs.  As seen below, the key/value pairs represent the name of the input on the left-hand side and the expression for the input's value on the right-hand side:

```
task hello {
  command {
    echo '${salutation} ${name}!'
  }
  output {
    String response = read_string(stdout())
  }
}

workflow test {
  call hello {
    input: salutation="greetings"
  }
  call hello as hello2
}
```

Now, the `hello.json` would require three inputs:

```
{
  "test.hello.name": "world",
  "test.hello2.name": "boston",
  "test.hello2.salutation": "hello"
}
```

Running this workflow still gives us the two greetings we expect:

```
$ java -jar cromwell.jar run hello.wdl hello.json
... truncated ...
{
  "test.hello.response": "greetings world!",
  "test.hello2.response": "hello boston!"
}
```

What if we wanted to parameterize the greeting and make it used for all invocations of task `hello`?  In this situation, a declaration can be used:

```
task hello {
  command {
    echo '${salutation}, ${name}!'
  }
  output {
    String response = read_string(stdout())
  }
}

workflow test {
  String greeting
  call hello {
    input: salutation=greeting
  }
  call hello as hello2 {
    input: salutation=greeting + " and nice to meet you"
  }
}
```

`String greeting` is referenced to satisfy the "salutation" parameter to both invocations of the `hello` task.

The inputs required to run this would be:

```
{
  "test.hello.name": "world",
  "test.hello2.name": "boston",
  "test.greeting": "hello"
}
```

And this would produce the following outputs when run

```
$ java -jar cromwell.jar run hello.wdl hello.json
... truncated ...
{
  "test.hello.response": "hello, world!",
  "test.hello2.response": "hello and nice to meet you, boston!"
}
```

## Using Files as Inputs

So far every example has used the default type of `String` for every input.  Passing files along to tasks is simply a matter of defining the input type as `File`:

```
task grep {
  command {
    grep -c '^...$' ${File file}
  }
  output {
    Int count = read_int(stdout())
  }
}

workflow test {
  call grep
}
```

The `read_int()` function here would read the contents of its parameter, and interpret the first line as an integer and return that value as a WDL `Int` type.

If I specified a file called `test_file` with the contents of:

```
foo
bar
baz
quux
```

And then the inputs JSON file would be:

```
{
  "test.grep.file": "test_file"
}
```

The result of running this would be:

```
$ java -jar cromwell.jar run grep.wdl grep.json
... truncated ...
{
  "test.grep.count": 3
}
```

## Scatter/Gather

Scatter blocks can be used to run the same call multiple times but only varying a specific parameter on each invocation.  Consider the following example:

```
task prepare {
  command <<<
    python -c "print('one\ntwo\nthree\nfour')"
  >>>
  output {
    Array[String] array = read_lines(stdout())
  }
}

task analysis {
  String str
  command <<<
    python -c "print('_${str}_')"
  >>>
  output {
    String out = read_string(stdout())
  }
}

task gather {
  Array[String] array
  command <<<
    echo ${sep=' ' array}
  >>>
  output {
    String str = read_string(stdout())
  }
}

workflow example {
  call prepare
  scatter (x in prepare.array) {
    call analysis {input: str=x}
  }
  call gather {input: array=analysis.out}
}
```

This example calls the `analysis` task once for each element in the array that the `prepare` task outputs.  The resulting outputs of this workflow would be:

```
{
  "example.analysis.out": ["_one_", "_two_", "_three_", "_four_"],
  "example.gather.str": "_one_ _two_ _three_ _four_",
  "example.prepare.array": ["one", "two", "three", "four"]
}
```

# Backends

A backend represents a way to run the user's command specified in the `task` section.  Currently three backends are supported:

* Local - Run jobs as subprocesses.  Supports launching in Docker containers.
* Sun GridEngine - Use `qsub` and job monitoring to run scripts.
* Google JES - Launch on Google Compute Cluster through JES.

Backends are specified via the configuration option `backend.backend` which can accept the values: sge, local, jes (e.g. `java -Dbackend.backend=sge`).

## Local Filesystem Assumptions / Layout

For backends `local` and `sge`, there are certain filesystem requirements that your environment needs to meet:

* The Cromwell process can write to the current working directory
* (`local` backend) Subprocesses that Cromwell launches can use child directories that Cromwell creates as their CWD.  The subprocess must have write access to the directory that Cromwell assigns as its current working directory.
* (`sge` backend) Jobs launched with `qsub` can use child directories as the working directory of the job, and write files to those directories.

The root directory that Cromwell uses for all workflows (`cromwell-root`) defaults to `./cromwell-executions`.  However, this is can be overwritten with the `-Dbackend.shared-filesystem.root=/your/path` option on the command line.

When cromwell runs a workflow through either of these two backends, it first creates a directory `<cromwell-root>/<workflow_uuid>`.  This is called the `workflow_root` and it is the root directory for all activity in this workflow.

Each `call` has its own subdirectory located at `<workflow_root>/call-<call_name>`.  This is the `<call_dir>`.  Within this directory are special files written by the backend and they're supposed to be backend specific things though there are commonalities.  For example, having a `stdout` and `stderr` file is common among both backends and they both write a shell script file to the `<call_dir>` as well.  See the descriptions below for details about backend-specific files that are written to these directories.

Any input files to the workflow need to be localized into `<workflow_root>/workflow-inputs`.  There are a few localization strategies that Cromwell will try until one works.  Below is the default order specified in `application.conf` but this order can be overridden:

* `hard-link` - This will create a hard link (not symbolic) link to the file
* `soft-link` - Create a symbolic link to the file.  This strategy is not applicable for tasks which specify a Docker image and will be ignored.
* `copy` - Make a copy the file

These options can be overridden with command line options to Java.  For instance, to use the strategies `copy` and `hard-link`, in that order:

```
java -Dbackend.shared-filesystem.localization.0=copy -Dbackend.shared-filesystem.localization.1=hard-link cromwell.jar ...
```

## Local Backend

The local backend will simply launch a subprocess and wait for it to exit.

This backend creates three files in the `<call_dir>` (see previous section):

* `script` - A shell script of the job to be run.  This contains the user's command from the `command` section of the WDL code.
* `stdout` - The standard output of the process
* `stderr` - The standard error of the process

The `script` file contains:

```
cd <container_call_root>
<user_command>
echo $? > rc
```

`<container_call_root>` would be equal to `<call_dir>` for non-Docker jobs, or it would be under `/root/<workflow_uuid>/call-<call_name>` if this is running in a Docker container.

The subprocess command that the local backend will launch is:

```
"/bin/bash" "-c" "cat script | <docker_run> /bin/bash <&0"
```

Where `<docker_run>` will be non-empty if this particular task specified a Docker container to run in.  `<docker_run>` looks like this:

```
docker run -v <local_workflow_dir>:/root/<workflow_uuid> -i <image>
```

## Sun GridEngine

The GridEngine backend uses `qsub` to launch a job and will poll the filesystem to determine if a job is completed.

This backend makes the same assumption about the filesystem that the local backend does: the Cromwell process and the jobs both have read/write access to the CWD of the job.

The CWD will contain a `script.sh` file which will contain:

```
\#!/bin/sh
<user_command>
echo $? > rc
```

The job is launched using the following command:

```
qsub -N <job_name> -V -b n -wd <call_dir> -o stdout -e stderr <call_dir>/script.sh
```

`<job_name>` is the string: `cromwell_<workflow_uuid_short>_<call_name>` (e.g. `cromwell_5103f8db_my_task`).

the `<call_dir>` contains the following special files added by the SGE backend:

* `qsub.stdout`, `qsub.stderr` - The results of the qsub command.
* `script.sh` - File containing the user's command and some wrapper code.
* `stdout`, `stderr` - Standard output streams of the actual job.
* `rc` - Return code of the SGE job, populated when the job has finished.

The SGE backend gets the job ID from parsing the `qsub.stdout` text file.

Since the `script.sh` ends with `echo $? > rc`, the backend will wait for the existence of this file, parse out the return code and determine success or failure and then subsequently post-process.

## Google JES

... more to come ...

# REST API

The `server` subcommand on the executable JAR will start an HTTP server which can accept WDL files to run as well as check status and output of existing workflows.

The following sub-sections define which HTTP Requests the web server can accept and what they will return.  Example HTTP requests are given in [HTTPie](https://github.com/jakubroztocil/httpie) and [cURL](http://curl.haxx.se/)

## REST API Versions

All web server requests include an API version in the url. The current version is `v1`.

## POST /workflows/:version

This endpoint accepts a POST request with a `multipart/form-data` encoded body.  The two elements in the body must be named `wdl` and `inputs`.  The `wdl` element contains the WDL file to run while the `inputs` contains a JSON file of the inputs to the workflow.

cURL:

```
$ curl -v "localhost:8000/workflows/v1" -F wdlSource=@src/main/resources/3step.wdl -F workflowInputs=@test.json
```

HTTPie:

```
$ http --print=hbHB --form POST localhost:8000/workflows/v1 wdlSource=@src/main/resources/3step.wdl workflowInputs@inputs.json
```

Request:

```
POST /workflows/v1 HTTP/1.1
Accept: */*
Accept-Encoding: gzip, deflate
Connection: keep-alive
Content-Length: 730
Content-Type: multipart/form-data; boundary=64128d499e9e4616adea7d281f695dca
Host: localhost:8000
User-Agent: HTTPie/0.9.2

--64128d499e9e4616adea7d281f695dca
Content-Disposition: form-data; name="wdlSource"

task ps {
  command {
    ps
  }
  output {
    File procs = stdout()
  }
}

task cgrep {
  command {
    grep '${pattern}' ${File in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
}

task wc {
  command {
    cat ${File in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
}

workflow three_step {
  call ps
  call cgrep {
    input: in_file=ps.procs
  }
  call wc {
    input: in_file=ps.procs
  }
}

--64128d499e9e4616adea7d281f695dca
Content-Disposition: form-data; name="workflowInputs"; filename="inputs.json"

{
    "three_step.cgrep.pattern": "..."
}

--64128d499e9e4616adea7d281f695dca--
```

Response:

```
HTTP/1.1 201 Created
Content-Length: 74
Content-Type: application/json; charset=UTF-8
Date: Tue, 02 Jun 2015 18:06:28 GMT
Server: spray-can/1.3.3

{
    "id": "69d1d92f-3895-4a7b-880a-82535e9a096e",
    "status": "Submitted"
}
```

## GET /workflows/:version/:id/status

cURL:

```
$ curl http://localhost:8000/workflows/v1/69d1d92f-3895-4a7b-880a-82535e9a096e/status
```

HTTPie:

```
$ http http://localhost:8000/workflows/v1/69d1d92f-3895-4a7b-880a-82535e9a096e/status
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 74
Content-Type: application/json; charset=UTF-8
Date: Tue, 02 Jun 2015 18:06:56 GMT
Server: spray-can/1.3.3

{
    "id": "69d1d92f-3895-4a7b-880a-82535e9a096e",
    "status": "Succeeded"
}
```

## GET /workflows/:version/:id/outputs

cURL:

```
$ curl http://localhost:8000/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/outputs
```

HTTPie:

```
$ http http://localhost:8000/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/outputs
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 241
Content-Type: application/json; charset=UTF-8
Date: Thu, 04 Jun 2015 12:15:33 GMT
Server: spray-can/1.3.3

{
    "id": "e442e52a-9de1-47f0-8b4f-e6e565008cf1",
    "outputs": {
        "three_step.cgrep.count": 8,
        "three_step.ps.procs": "/var/folders/kg/c7vgxnn902lc3qvc2z2g81s89xhzdz/T/stdout2814345504446060277.tmp",
        "three_step.wc.count": 8
    }
}
```


## GET /workflows/:version/:id/outputs/:call

cURL:

```
$ curl http://localhost:8000/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/outputs/three_step.wc
```

HTTPie:

```
$ http http://localhost:8000/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/outputs/three_step.wc
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 241
Content-Type: application/json; charset=UTF-8
Date: Thu, 04 Jun 2015 12:15:33 GMT
Server: spray-can/1.3.3

{
    "id": "e442e52a-9de1-47f0-8b4f-e6e565008cf1",
    "outputs": {
        "three_step.wc.count": 8
    }
}
```
## GET /workflows/:version/:id/logs/:call

This will return paths to the standard out and standard error files that were generated during the execution
of a particular fully-qualified name for a call.

cURL:

```
$ curl http://localhost:8000/workflows/v1/b3e45584-9450-4e73-9523-fc3ccf749848/logs/three_step.wc
```

HTTPie:

```
$ http http://localhost:8000/workflows/v1/b3e45584-9450-4e73-9523-fc3ccf749848/logs/three_step.wc
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 379
Content-Type: application/json; charset=UTF-8
Date: Mon, 03 Aug 2015 17:11:28 GMT
Server: spray-can/1.3.3

{
    "id": "b3e45584-9450-4e73-9523-fc3ccf749848",
    "logs": {
        "three_step.wc": {
            "stderr": "/Users/sfrazer/projects/cromwell/cromwell-executions/test/b3e45584-9450-4e73-9523-fc3ccf749848/three_step.wc/stderr6126967977036995110.tmp",
            "stdout": "/Users/sfrazer/projects/cromwell/cromwell-executions/test/b3e45584-9450-4e73-9523-fc3ccf749848/three_step.wc/stdout6128485235785447571.tmp"
        }
    }
}
```
## GET /workflows/:version/:id/logs

cURL:

```
$ curl http://localhost:8000/workflows/v1/b3e45584-9450-4e73-9523-fc3ccf749848/logs
```

HTTPie:

```
$ http http://localhost:8000/workflows/v1/b3e45584-9450-4e73-9523-fc3ccf749848/logs
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 379
Content-Type: application/json; charset=UTF-8
Date: Mon, 03 Aug 2015 17:11:28 GMT
Server: spray-can/1.3.3

{
    "id": "b3e45584-9450-4e73-9523-fc3ccf749848",
    "logs": {
        "call.ps": {
            "stderr": "/Users/sfrazer/projects/cromwell/cromwell-executions/test/b3e45584-9450-4e73-9523-fc3ccf749848/call-ps/stderr6126967977036995110.tmp",
            "stdout": "/Users/sfrazer/projects/cromwell/cromwell-executions/test/b3e45584-9450-4e73-9523-fc3ccf749848/call-ps/stdout6128485235785447571.tmp"
        },
        "call.cgrep": {
            "stderr": "/Users/sfrazer/projects/cromwell/cromwell-executions/test/b3e45584-9450-4e73-9523-fc3ccf749848/call-cgrep/stderr6126967977036995110.tmp",
            "stdout": "/Users/sfrazer/projects/cromwell/cromwell-executions/test/b3e45584-9450-4e73-9523-fc3ccf749848/call-cgrep/stdout6128485235785447571.tmp"
        },
        "call.wc": {
            "stderr": "/Users/sfrazer/projects/cromwell/cromwell-executions/test/b3e45584-9450-4e73-9523-fc3ccf749848/call-wc/stderr6126967977036995110.tmp",
            "stdout": "/Users/sfrazer/projects/cromwell/cromwell-executions/test/b3e45584-9450-4e73-9523-fc3ccf749848/call-wc/stdout6128485235785447571.tmp"
        }
    }
}
```

## POST /workflows/:version/:id/abort

cURL:

```
$ curl -X POST http://localhost:8000/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/abort
```

HTTPie:

```
$ http POST http://localhost:8000/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/abort
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 241
Content-Type: application/json; charset=UTF-8
Date: Thu, 04 Jun 2015 12:15:33 GMT
Server: spray-can/1.3.3

{
    "id": "e442e52a-9de1-47f0-8b4f-e6e565008cf1",
    "status": "Aborted"
}
```

# Developer

## Generate WDL Parser

Install the latest version of [Hermes](http://github.com/scottfrazer/hermes), then run the following command within this directory:

```
hermes generate src/main/resources/grammar.hgr --language=java --directory=src/main/java --name=wdl --java-package=cromwell.parser --java-use-apache-commons --header
```

The grammar for the WDL lexer/parser is defined in `src/main/resources/grammar.hgr`.  Any changes to that grammar should result in a regeneration of the parser and then run the unit tests.  Changing the AST could be disruptive if keys are renamed or objects restructured too much.  It's best to find these issues as soon as possible.

## Generating and Hosting ScalaDoc

Essentially run `sbt doc` then commit the generated code into the `gh-pages` branch on this repository

```
$ sbt doc
$ git co gh-pages
$ mv target/scala-2.11/api scaladoc
$ git add scaladoc
$ git commit -m "API Docs"
$ git push origin gh-pages
```

## Architecture

![Cromwell Architecture](http://i.imgur.com/kPPTe0l.png)

The architecture is split into four layers, from bottom to top:

### cromwell.parser

Contains only the WDL parser to convert WDL source code to an abstract syntax tree.  Clients should never need to interact with WDL at this level, though nothing specifically precludes that.

### cromwell.binding

Contains code that takes an abstract syntax tree and returns native Scala object representations of those ASTs.  This layer will also have functions for evaluating expressions when support for that is added.

### cromwell.backend

Contains implementations of an interface to launch jobs.  `cromwell.engine` will use this to execute and monitor jobs.

### cromwell.engine

![Engine Actors](http://i.imgur.com/sF9vMt2.png)

Contains the Akka code and actor system to execute a workflow.  This layer should operate entirely on objects returned from the `cromwell.binding` layer.
