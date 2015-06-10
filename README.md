[![Build Status](https://travis-ci.org/broadinstitute/cromwell.svg?branch=develop)](https://travis-ci.org/broadinstitute/cromwell?branch=develop)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/cromwell/badge.svg?branch=develop)](https://coveralls.io/r/broadinstitute/cromwell?branch=develop)


# Cromwell

Workflow engine using [WDL](https://github.com/broadinstitute/wdl/blob/wdl2/SPEC.md) as the workflow and task language.

<!---toc start-->

* [Cromwell](#cromwell)
  * [API Documenation](#api-documenation)
  * [Scala API Usage](#scala-api-usage)
  * [REST API](#rest-api)
    * [POST /workflows](#post-workflows)
    * [GET /workflow/:id/status](#get-workflowidstatus)
    * [GET /workflow/:id/outputs](#get-workflowidoutputs)
* [Developer](#developer)
  * [Generate WDL Parser](#generate-wdl-parser)
  * [Generating and Hosting ScalaDoc](#generating-and-hosting-scaladoc)
  * [Architecture](#architecture)
    * [cromwell.parser](#cromwellparser)
    * [cromwell.binding](#cromwellbinding)
    * [cromwell.backend](#cromwellbackend)
    * [cromwell.engine](#cromwellengine)

<!---toc end-->

## API Documenation

API Documentation can be found [here](http://broadinstitute.github.io/cromwell/scaladoc).

Note that this may not be completely up to date or even useful at this time.

## Scala API Usage

The main entry point into the parser is the `WdlNamespace` object.  A WDL file is considered a namespace, and other namespaces can be included by using the `import` statement (but only with an `as` clause).

```scala
import cromwell.binding.WdlNamespace

object Main extends App {
  val ns = WdlNamespace.load(new File("test.wdl"))
  val ns2 = WdlNamespace.load("workflow wf {}")
  val workflows = ns.workflows
  val tasks = ns.tasks
  val commands = ns.tasks.map{_.command}
}
```

To access only the parser, use the `AstTools` library, as follows:

```scala
import cromwell.binding.WdlNamespace

object Main extends App {
  val ast = AstTools.getAst(new File(args(0)))

  // Second parameter is a descriptor about where the first string came from.
  // Most of the time this would be the URI of where the text was loaded from,
  // but there are no restrictions on what the string can be.
  val astFromString = AstTools.getAst("workflow wf {}", "string")

  println(ast.toPrettyString)
}
```

## Command Line Usage

Run the JAR file with no arguments to get the usage message:

```
$ java -jar target/scala-2.11/cromwell-0.1-SNAPSHOT.jar

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

server

  Starts a web server on port 8000.  See the web server
  documentation for more details about the API endpoints.
```

### validate

Given a WDL file, this runs the full syntax checker over the file and resolves imports in the process.  If any syntax errors are found, they are printed out.  Otherwise the program exits.

Error if a `call` references a task that doesn't exist:

```
$ java -jar ../target/scala-2.11/cromwell-0.1-SNAPSHOT.jar validate 2.wdl
ERROR: Call references a task (BADps) that doesn't exist (line 22, col 8)

  call BADps
       ^

```

Error if namespace and task have the same name:

```
$ java -jar ../target/scala-2.11/cromwell-0.1-SNAPSHOT.jar validate 5.wdl
ERROR: Task and namespace have the same name:

Task defined here (line 3, col 6):

task ps {
     ^

Import statement defined here (line 1, col 20):

import "ps.wdl" as ps
                   ^
```

### inputs

Examine a WDL file with one workflow in it, compute all the inputs needed for that workflow and output a JSON template that the user can fill in with values.  The keys in this document should remain unchanged.  The values tell you what type the parameter is expecting.  For example, if the value were `Array[String]`, then it's expecting a JSON array of JSON strings, like this: `["string1", "string2", "string3"]`

```
$ java -jar ../target/scala-2.11/cromwell-0.1-SNAPSHOT.jar inputs 3step.wdl
{
  "three_step.cgrep.pattern": "String"
}
```

This inputs document is used as input to the `run` subcommand.

### run

Given a WDL file and a JSON inputs file (see `inputs` subcommand), Run the workflow and print the outputs:

```
$ java -jar ../target/scala-2.11/cromwell-0.1-SNAPSHOT.jar run 3step.wdl inputs.json
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

### parse

Given a WDL file input, this does grammar level syntax checks and prints out the resulting abstract syntax tree.

```
$ echo "workflow wf {}" | java -jar ../target/scala-2.11/cromwell-0.1-SNAPSHOT.jar parse /dev/stdin
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

### server

Start a server on port 8000, the REST API for this server is described in the next section

## REST API

The `server` subcommand on the executable JAR will start an HTTP server which can accept WDL files to run as well as check status and output of existing workflows.

The following sub-sections define which HTTP Requests the web server can accept and what they will return.

### POST /workflows

This endpoint accepts a POST request with a `multipart/form-data` encoded body.  The two elements in the body must be named `wdl` and `inputs`.  The `wdl` element contains the WDL file to run while the `inputs` contains a JSON file of the inputs to the workflow.

```
$ http --print=hbHB --form POST localhost:8000/workflows wdlSource=@src/main/resources/3step.wdl workflowInputs@inputs.json
```

Request:

```
POST /workflows HTTP/1.1
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

### GET /workflow/:id/status

```
http http://localhost:8000/workflow/69d1d92f-3895-4a7b-880a-82535e9a096e/status
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

### GET /workflow/:id/outputs

```
http http://localhost:8000/workflow/e442e52a-9de1-47f0-8b4f-e6e565008cf1/outputs
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

# Developer

## Generate WDL Parser

Install the latest version of [Hermes](http://github.com/scottfrazer/hermes), then run the following command within this directory:

```
hermes generate src/main/resources/grammar.hgr --language=java --directory=src/main/java --name=wdl --java-package=cromwell.parser --java-use-apache-commons
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

![Engine Actors](http://i.imgur.com/ByFUakW.png)

Contains the Akka code and actor system to execute a workflow.  This layer should operate entirely on objects returned from the `cromwell.binding` layer.

|Start|End|Message|Parameters|State Change|
|-----|---|-------|----------|------------|
|![WMA](http://i.imgur.com/98ugXkZ.png)|![WA](http://i.imgur.com/eJxn3wu.png)|Start||
|![WMA](http://i.imgur.com/98ugXkZ.png)|![WA](http://i.imgur.com/eJxn3wu.png)|SubscribeTransitionCallback|WorkflowManagerActor|
|![WA](http://i.imgur.com/eJxn3wu.png)|![SA](http://i.imgur.com/JzNNe64.png)|UpdateStatus|Call,status|
|![WA](http://i.imgur.com/eJxn3wu.png)|![SA](http://i.imgur.com/JzNNe64.png)|CallCompleted|Call,outputs|
|![WA](http://i.imgur.com/eJxn3wu.png)|![SA](http://i.imgur.com/JzNNe64.png)|GetOutputs||
|![WA](http://i.imgur.com/eJxn3wu.png)|![CA](http://i.imgur.com/Nyoln74.png)|Start||

> **TODO**: This table is not complete.

