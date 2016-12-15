[![Build Status](https://travis-ci.org/broadinstitute/cromwell.svg?branch=develop)](https://travis-ci.org/broadinstitute/cromwell?branch=develop)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/cromwell/badge.svg?branch=develop)](https://coveralls.io/r/broadinstitute/cromwell?branch=develop)
[![Join the chat at https://gitter.im/broadinstitute/cromwell](https://badges.gitter.im/broadinstitute/cromwell.svg)](https://gitter.im/broadinstitute/cromwell?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

Cromwell
========

A [Workflow Management System](https://en.wikipedia.org/wiki/Workflow_management_system) geared towards [scientific workflows](https://en.wikipedia.org/wiki/Scientific_workflow_system). Cromwell is open sourced under the BSD 3-Clause license.

<!---toc start-->

* [Getting Help](#getting-help)
  * [Website and User Guide](#website-and-user-guide)
  * [Support Forum](#support-forum)
  * [Gitter](#gitter)
* [Requirements](#requirements)
* [Building](#building)
* [Installing](#installing)
  * [Upgrading from 0.19 to 0.21](#upgrading-from-019-to-021)
* [Command Line Usage](#command-line-usage)
  * [run](#run)
  * [server](#server)
  * [version](#version)
* [Getting Started with WDL](#getting-started-with-wdl)
  * [WDL Support](#wdl-support)
* [Configuring Cromwell](#configuring-cromwell)
  * [Workflow Submission](#workflow-submission)
  * [Database](#database)
  * [SIGINT abort handler](#sigint-abort-handler)
* [Security](#security)
* [Backends](#backends)
  * [Backend Filesystems](#backend-filesystems)
    * [Shared Local Filesystem](#shared-local-filesystem)
    * [Google Cloud Storage Filesystem](#google-cloud-storage-filesystem)
  * [Local Backend](#local-backend)
  * [Sun GridEngine Backend](#sun-gridengine-backend)
  * [HtCondor Backend](#htcondor-backend)
    * [Caching configuration](#caching-configuration)
    * [Docker](#docker)
    * [CPU, Memory and Disk](#cpu-memory-and-disk)
    * [Native Specifications](#native-specifications)
  * [Spark Backend](#spark-backend)
    * [Configuring Spark Project](#configuring-spark-project)
    * [Configuring Spark Master and Deploy Mode](#configuring-spark-master-and-deploy-mode)
    * [Spark runtime attributes](#spark-runtime-attributes)
    * [Spark Environment](#spark-environment)
    * [Sample Wdl](#sample-wdl)
  * [Google JES Backend](#google-jes-backend)
    * [Configuring Google Project](#configuring-google-project)
    * [Configuring Authentication](#configuring-authentication)
      * [Application Default Credentials](#application-default-credentials)
      * [Service Account](#service-account)
      * [Refresh Token](#refresh-token)
    * [Docker](#docker)
    * [Monitoring](#monitoring)
* [Runtime Attributes](#runtime-attributes)
  * [Specifying Default Values](#specifying-default-values)
  * [continueOnReturnCode](#continueonreturncode)
  * [cpu](#cpu)
  * [disks](#disks)
    * [Boot Disk](#boot-disk)
  * [zones](#zones)
  * [docker](#docker)
  * [failOnStderr](#failonstderr)
  * [memory](#memory)
  * [preemptible](#preemptible)
* [Logging](#logging)
* [Workflow Options](#workflow-options)
* [Call Caching](#call-caching)
  * [Configuring Call Caching](#configuring-call-caching)
  * [Call Caching Workflow Options](#call-caching-workflow-options)
  * [Local Filesystem Options](#local-filesystem-options)
* [Imports](#imports)
* [Sub Workflows](#sub-workflows)
* [Meta blocks](#meta-blocks)
* [REST API](#rest-api)
  * [REST API Versions](#rest-api-versions)
  * [POST /api/workflows/:version](#post-apiworkflowsversion)
  * [POST /api/workflows/:version/batch](#post-apiworkflowsversionbatch)
  * [GET /api/workflows/:version/query](#get-apiworkflowsversionquery)
  * [POST /api/workflows/:version/query](#post-apiworkflowsversionquery)
  * [GET /api/workflows/:version/:id/status](#get-apiworkflowsversionidstatus)
  * [GET /api/workflows/:version/:id/outputs](#get-apiworkflowsversionidoutputs)
  * [GET /api/workflows/:version/:id/timing](#get-apiworkflowsversionidtiming)
  * [GET /api/workflows/:version/:id/logs](#get-apiworkflowsversionidlogs)
  * [GET /api/workflows/:version/:id/metadata](#get-apiworkflowsversionidmetadata)
  * [POST /api/workflows/:version/:id/abort](#post-apiworkflowsversionidabort)
  * [GET /api/workflows/:version/backends](#get-apiworkflowsversionbackends)
  * [GET /api/engine/:version/stats](#get-apiengineversionstats)
  * [GET /api/engine/:version/version](#get-apiengineversionversion)
  * [Error handling](#error-handling)
* [Developer](#developer)
  * [Generating table of contents on Markdown files](#generating-table-of-contents-on-markdown-files)
  * [Generating and Hosting ScalaDoc](#generating-and-hosting-scaladoc)

<!---toc end-->

# Getting Help

## Website and User Guide

The [WDL website](https://software.broadinstitute.org/wdl/) is the best place to go for more information on both WDL and Cromwell. In particular new users should check out the [user guide](https://software.broadinstitute.org/wdl/userguide/) which has many tutorials, examples and other bits to get you started.

## Support Forum

If you have questions that aren't covered by the website you can ask them in the [support forum](http://gatkforums.broadinstitute.org/wdl/categories/ask-the-wdl-team).

## Gitter
There is a [Cromwell gitter channel](https://gitter.im/broadinstitute/cromwell) where people can discuss Cromwell and related topics with both the developers and user community.

# Requirements

The following is the toolchain used for development of Cromwell.  Other versions may work, but these are recommended.

* [Scala 2.11.8](http://www.scala-lang.org/news/2.11.8/)
* [SBT 0.13.12](https://github.com/sbt/sbt/releases/tag/v0.13.12)
* [Java 8](http://www.oracle.com/technetwork/java/javase/overview/java8-2100321.html)

# Building

`sbt assembly` will build a runnable JAR in `target/scala-2.11/`

Tests are run via `sbt test`.  Note that the tests do require Docker to be running.  To test this out while downloading the Ubuntu image that is required for tests, run `docker pull ubuntu:latest` prior to running `sbt test`

# Installing

OS X users can install Cromwell with Homebrew: `brew install cromwell`.

## Upgrading from 0.19 to 0.21

See the [migration document](MIGRATION.md) for more details.

# Command Line Usage

Run the JAR file with no arguments to get the usage message:

```


$ java -jar cromwell.jar
java -jar cromwell.jar <action> <parameters>

Actions:
run <WDL file> [<JSON inputs file>] [<JSON workflow options>]
  [<OUTPUT workflow metadata>] [<Zip of WDL Files>]

  Given a WDL file and JSON file containing the value of the
  workflow inputs, this will run the workflow locally and
  print out the outputs in JSON format.  The workflow
  options file specifies some runtime configuration for the
  workflow (see README for details).  The workflow metadata
  output is an optional file path to output the metadata. The
  directory of WDL files is optional. However, it is required
  if the primary workflow imports workflows that are outside
  of the root directory of the Cromwell project.

  Use a single dash ("-") to skip optional files. Ex:
    run noinputs.wdl - - metadata.json -

  server

  Starts a web server on port 8000.  See the web server
  documentation for more details about the API endpoints.

  -version

  Returns the version of the Cromwell engine.

```

## run

Given a WDL file and a JSON inputs file (see `inputs` subcommand), Run the workflow and print the outputs:

```
$ java -jar cromwell.jar run 3step.wdl inputs.json
... play-by-play output ...
{
  "three_step.ps.procs": "/var/folders/kg/c7vgxnn902lc3qvc2z2g81s89xhzdz/T/stdout1272284837004786003.tmp",
  "three_step.cgrep.count": 0,
  "three_step.wc.count": 13
}
```

The JSON inputs can be left off if there's a file with the same name as the WDL file but with a `.inputs` extension.  For example, this will assume that `3step.inputs` exists:

```
$ java -jar cromwell.jar run 3step.wdl
```

If your workflow has no inputs, you can specify `-` as the value for the inputs parameter:

```
$ java -jar cromwell.jar run my_workflow.wdl -
```

The third, optional parameter to the 'run' subcommand is a JSON file of workflow options.  By default, the command line will look for a file with the same name as the WDL file but with the extension `.options`.  But one can also specify a value of `-` manually to specify that there are no workflow options.

See the section [workflow options](#workflow-options) for more details.

```
$ java -jar cromwell.jar run my_jes_wf.wdl my_jes_wf.json wf_options.json
```

The fourth, optional parameter to the 'run' subcommand is a path where the workflow metadata will be written.  By default, no workflow metadata will be written.

```
$ java -jar cromwell.jar run my_wf.wdl - - my_wf.metadata.json
... play-by-play output ...
$ cat my_wf.metadata.json
{
  "workflowName": "w",
  "calls": {
    "w.x": [{
      "executionStatus": "Done",
      "stdout": "/Users/jdoe/projects/cromwell/cromwell-executions/w/a349534f-137b-4809-9425-1893ac272084/call-x/stdout",
      "shardIndex": -1,
      "outputs": {
        "o": "local\nremote"
      },
      "runtimeAttributes": {
        "failOnStderr": "false",
        "continueOnReturnCode": "0"
      },
      "cache": {
        "allowResultReuse": true
      },
      "inputs": {
        "remote": "/Users/jdoe/remote.txt",
        "local": "local.txt"
      },
      "returnCode": 0,
      "backend": "Local",
      "end": "2016-07-11T10:27:56.074-04:00",
      "stderr": "/Users/jdoe/projects/cromwell/cromwell-executions/w/a349534f-137b-4809-9425-1893ac272084/call-x/stderr",
      "callRoot": "cromwell-executions/w/a349534f-137b-4809-9425-1893ac272084/call-x",
      "attempt": 1,
      "start": "2016-07-11T10:27:55.992-04:00"
    }]
  },
  "outputs": {
    "w.x.o": "local\nremote"
  },
  "workflowRoot": "cromwell-executions/w/a349534f-137b-4809-9425-1893ac272084",
  "id": "a349534f-137b-4809-9425-1893ac272084",
  "inputs": {
    "w.x.remote": "/Users/jdoe/remote.txt",
    "w.x.local": "local.txt"
  },
  "submission": "2016-07-11T10:27:54.907-04:00",
  "status": "Succeeded",
  "end": "2016-07-11T10:27:56.108-04:00",
  "start": "2016-07-11T10:27:54.919-04:00"
}
```

The fifth, optional parameter to the 'run' subcommand is a zip file which contains WDL source files. This zip file can be passed
and your primary workflow can import any WDL's from that collection and re-use those tasks.

For example, consider you have a directory of WDL files:
```
my_WDLs
└──cgrep.wdl
└──ps.wdl
└──wc.wdl
```

If you zip that directory to my_WDLs.zip, you have the option to pass it in as the last parameter in your run command
and be able to reference these WDLs as imports in your primary WDL. For example, your primary WDL can look like this:
```
import "ps.wdl" as ps
import "cgrep.wdl"
import "wc.wdl" as wordCount

workflow threestep {

call ps.ps as getStatus
call cgrep.cgrep { input: str = getStatus.x }
call wordCount { input: str = ... }

}

```
The command to run this WDL, without needing any inputs, workflow options or metadata files would look like:

```
$ java -jar cromwell.jar run threestep.wdl - - - /path/to/my_WDLs.zip
```

## server

Start a server on port 8000, the API for the server is described in the [REST API](#rest-api) section.

## version

Returns the version of Cromwell engine.

# Getting Started with WDL

For many examples on how to use WDL see [the WDL site](https://github.com/broadinstitute/wdl#getting-started-with-wdl)

## WDL Support

:pig2: Cromwell supports the following subset of WDL features:

* [Language Specification](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#language-specification)
  * [Whitespace, Strings, Identifiers, Constants](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#whitespace-strings-identifiers-constants)
  * [Types](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#types)
  * [Fully Qualified Names & Namespaced Identifiers](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#fully-qualified-names--namespaced-identifiers)
  * [Declarations](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#declarations)
  * [Expressions](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#expressions)
  * [Operator Precedence Table](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#operator-precedence-table)
  * [Member Access](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#member-access)
  * [Map and Array Indexing](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#map-and-array-indexing)
  * [Function Calls](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#function-calls)
  * [Array Literals](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#array-literals)
  * [Map Literals](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#map-literals)
* [Document](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#document)
* [Task Definition](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#task-definition)
  * [Sections](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#sections)
  * [Command Section](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#command-section)
    * [Command Parts](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#command-parts)
    * [Command Part Options](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#command-part-options)
      * [sep](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#sep)
      * [true and false](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#true-and-false)
      * [default](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#default)
    * [Alternative heredoc syntax](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#alternative-heredoc-syntax)
    * [Stripping Leading Whitespace](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#stripping-leading-whitespace)
  * [Outputs Section](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#outputs-section)
  * [String Interpolation](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#string-interpolation)
  * [Runtime Section](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#runtime-section)
    * [docker](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#docker)
    * [memory](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#memory)
  * [Parameter Metadata Section](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#parameter-metadata-section)
  * [Metadata Section](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#metadata-section)
* [Workflow Definition](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#workflow-definition)
  * [Call Statement](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#call-statement)
    * [Sub Workflows](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#sub-workflows)
  * [Scatter](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#scatter)
  * [Parameter Metadata](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#parameter-metadata)
  * [Metadata](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#metadata)
  * [Outputs](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#outputs)
* [Namespaces](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#namespaces)
* [Scope](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#scope)
* [Optional Parameters & Type Constraints](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#optional-parameters--type-constraints)
  * [Prepending a String to an Optional Parameter](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#prepending-a-string-to-an-optional-parameter)
* [Scatter / Gather](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#scatter--gather)
* [Variable Resolution](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#variable-resolution)
  * [Task-Level Resolution](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#task-level-resolution)
  * [Workflow-Level Resolution](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#workflow-level-resolution)
* [Computing Inputs](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#computing-inputs)
  * [Task Inputs](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#task-inputs)
  * [Workflow Inputs](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#workflow-inputs)
  * [Specifying Workflow Inputs in JSON](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#specifying-workflow-inputs-in-json)
* [Type Coercion](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#type-coercion)
* [Standard Library](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#standard-library)
  * [File stdout()](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#file-stdout)
  * [File stderr()](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#file-stderr)
  * [Array\[String\] read_lines(String|File)](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#arraystring-read_linesstringfile)
  * [Array\[Array\[String\]\] read_tsv(String|File)](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#arrayarraystring-read_tsvstringfile)
  * [Map\[String, String\] read_map(String|File)](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#mapstring-string-read_mapstringfile)
  * [Object read_object(String|File)](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#object-read_objectstringfile)
  * [Array\[Object\] read_objects(String|File)](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#arrayobject-read_objectsstringfile)
  * [Int read_int(String|File)](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#int-read_intstringfile)
  * [String read_string(String|File)](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#string-read_stringstringfile)
  * [Float read_float(String|File)](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#float-read_floatstringfile)
  * [Boolean read_boolean(String|File)](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#boolean-read_booleanstringfile)
  * [File write_lines(Array\[String\])](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#file-write_linesarraystring)
  * [File write_tsv(Array\[Array\[String\]\])](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#file-write_tsvarrayarraystring)
  * [File write_map(Map\[String, String\])](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#file-write_mapmapstring-string)
  * [File write_object(Object)](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#file-write_objectobject)
  * [File write_objects(Array\[Object\])](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#file-write_objectsarrayobject)
  * [Float size(File, \[String\])](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#float-sizefile-string)
  * [String sub(String, String, String)](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#string-substring-string-string)
  * [Array\[Int\] range(Int)](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#arrayint-rangeint)
  * [Array\[Array\[X\]\] transpose(Array\[Array\[X\]\])](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#arrayarrayx-transposearrayarrayx)
  * [Pair(X,Y) zip(X,Y)](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#pairxy-zipxy)
  * [Pair(X,Y) cross(X,Y)](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#pairxy-crossxy)
* [Data Types & Serialization](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#data-types--serialization)
  * [Serialization of Task Inputs](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#serialization-of-task-inputs)
    * [Primitive Types](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#primitive-types)
    * [Compound Types](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#compound-types)
      * [Array serialization](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#array-serialization)
        * [Array serialization by expansion](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#array-serialization-by-expansion)
        * [Array serialization using write_lines()](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#array-serialization-using-write_lines)
      * [Map serialization](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#map-serialization)
        * [Map serialization using write_map()](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#map-serialization-using-write_map)
      * [Object serialization](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#object-serialization)
        * [Object serialization using write_object()](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#object-serialization-using-write_object)
      * [Array\[Object\] serialization](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#arrayobject-serialization)
        * [Array\[Object\] serialization using write_objects()](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#arrayobject-serialization-using-write_objects)
  * [De-serialization of Task Outputs](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#de-serialization-of-task-outputs)
    * [Primitive Types](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#primitive-types)
    * [Compound Types](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#compound-types)
      * [Array deserialization](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#array-deserialization)
        * [Array deserialization using read_lines()](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#array-deserialization-using-read_lines)
      * [Map deserialization](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#map-deserialization)
        * [Map deserialization using read_map()](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#map-deserialization-using-read_map)
      * [Object deserialization](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#object-deserialization)
        * [Object deserialization using read_object()](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#object-deserialization-using-read_object)
      * [Array\[Object\] deserialization](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#arrayobject-deserialization)
        * [Object deserialization using read_objects()](https://github.com/broadinstitute/wdl/blob/develop/SPEC.md#object-deserialization-using-read_objects)


# Configuring Cromwell

Cromwell's default configuration file is located at `src/main/resources/application.conf`.

The configuration file is in [Hocon](https://github.com/typesafehub/config/blob/master/HOCON.md#hocon-human-optimized-config-object-notation) which means the configuration file can specify configuration as JSON-like stanzas like:

```hocon
webservice {
  port = 8000
  interface = 0.0.0.0
  instance.name = "reference"
}
```

Or, alternatively, as dot-separated values:

```hocon
webservice.port = 8000
webservice.interface = 0.0.0.0
webservice.instance.name = "reference"
```

This allows any value to be overridden on the command line:

```
java -Dwebservice.port=8080 cromwell.jar ...
```

It is recommended that one copies `src/main/resources/application.conf`, modify it, then link to it via:

```
java -Dconfig.file=/path/to/application.conf cromwell.jar ...
```

## Workflow Submission

Cromwell has a configurable cap on the number of workflows running at a time. To set this value provide an integer value to the `system.max-concurrent-workflows` config value.

Cromwell will look for new workflows to start on a regular interval which can be modified by setting the `system.new-workflow-poll-rate` config value, which is the number of seconds between workflow launches. On every poll, Cromwell will take at most `system.max-workflow-launch-count` new submissions, provided there are new workflows to launch and the `system.max-concurrent-workflows` number has not been reached.

## Database

Cromwell uses either an in-memory or MySQL database to track the execution of workflows and store outputs of task invocations.

By default, Cromwell uses an in-memory database which will only live for the duration of the JVM.  This provides a quick way to run workflows locally without having to set up MySQL, though it also makes workflow executions somewhat transient.

To configure Cromwell to instead point to a MySQL database, first create the empty database.  In the example below, the database name is `cromwell`.

Then, edit the configuration file `database` stanza, as follows:

```
database {

  driver = "slick.driver.MySQLDriver$"
  db {
    driver = "com.mysql.jdbc.Driver"
    url = "jdbc:mysql://host/cromwell"
    user = "user"
    password = "pass"
    connectionTimeout = 5000
  }

  test {
    ...
  }
}
```

## SIGINT abort handler

For backends that support aborting task invocations, Cromwell can be configured to automatically try to abort all currently running calls (and set their status to `Aborted`) when a SIGINT is sent to the Cromwell process.  To turn this feature on, set the configuration option

```
system {
  abort-jobs-on-terminate=true
}
```

Or, via `-Dsystem.abort-jobs-on-terminate=true` command line option.

By default, this value is false when running `java -jar cromwell.jar server`, and true when running `java -jar cromwell.jar run <wdl> <inputs>`.

# Security

 - Cromwell is NOT on its own a security appliance!
 - Only YOU are responsible for your own security! 
 - Some recommendations and suggestions on security can be found in the [SecurityRecommendations.md](SecurityRecommendations.md) document

# Backends

A backend represents a way to run the user's command specified in the `task` section.  Cromwell allows for backends conforming to
the Cromwell backend specification to be plugged into the Cromwell engine.  Additionally, backends are included with the
Cromwell distribution:

* Local / GridEngine / LSF / etc. - Run jobs as subprocesses or via a dispatcher.  Supports launching in Docker containers. Use `bash`, `qsub`, `bsub`, etc. to run scripts.
* Google JES - Launch jobs on Google Compute Engine through the Job Execution Service (JES).
* HtCondor - Allows to execute jobs using HTCondor.
* Spark - Adds support for execution of spark jobs.

Backends are specified in the `backend` configuration block under `providers`.  Each backend has a configuration that looks like:

```hocon
backend {
  default = "Local"
  providers {
    BackendName {
      actor-factory = "FQN of BackendLifecycleActorFactory instance"
      config {
        key = "value"
        key2 = "value2"
        ...
      }
    }
  }
}
```

The structure within the `config` block will vary from one backend to another; it is the backend implementation's responsibility
to be able to interpret its configuration.

In the example below two backend types are named within the `providers` section here, so both
are available.  The default backend is specified by `backend.default` and must match the `name` of one of the
configured backends:

```hocon
backend {
  default = "Local"
  providers {
    Local {
      actor-factory = "cromwell.backend.impl.local.LocalBackendLifecycleActorFactory"
      config {
        root: "cromwell-executions"
        filesystems = {
          local {
            localization: [
              "hard-link", "soft-link", "copy"
            ]
          }
          gcs {
            # References an auth scheme defined in the 'google' stanza.
            auth = "application-default"
          }
        }
      }
    },
    JES {
      actor-factory = "cromwell.backend.impl.jes.JesBackendLifecycleActorFactory"
      config {
        project = "my-cromwell-workflows"
        root = "gs://my-cromwell-workflows-bucket"
        maximum-polling-interval = 600
        dockerhub {
          # account = ""
          # token = ""
        }
        genomics {
          # A reference to an auth defined in the 'google' stanza at the top.  This auth is used to create
          # Pipelines and manipulate auth JSONs.
          auth = "application-default"
          endpoint-url = "https://genomics.googleapis.com/"
        }
        filesystems = {
          gcs {
            # A reference to a potentially different auth for manipulating files via engine functions.
            auth = "user-via-refresh"
          }
        }
      }
    }
  ]
}
```

## Backend Filesystems

Each backend will utilize filesystems to store the directory structure of an executed workflow.  Currently, the backends and the type of filesystems that the backend use are tightly coupled.  In future versions of Cromwell, they may be more loosely coupled.

The backend/filesystem pairings are as follows:

* [Local Backend](#local-backend) and associated backends primarily use the [Shared Local Filesystem](#shared-local-filesystem).
* [JES Backend](#google-jes-backend) uses the [Google Cloud Storage Filesystem](#google-cloud-storage-filesystem)

Note that while Local, SGE, LSF, etc. backends use the local or network filesystem for the directory structure of a workflow, they are able to localize inputs
from GCS paths if configured to use a GCS filesystem.  See [Google Cloud Storage Filesystem](#google-cloud-storage-filesystem) for more details.

### Shared Local Filesystem

For the [local](#local-backend) and [Sun GridEngine](#sun-gridengine-backend) backends, the following is required of the underlying filesystem:

Cromwell is configured with a root execution directory which is set in the configuration file under `backend.providers.<backend_name>.config.root`.  This is called the `cromwell_root` and it is set to `./cromwell-executions` by default.  Relative paths are interpreted as relative to the current working directory of the Cromwell process.

When Cromwell runs a workflow, it first creates a directory `<cromwell_root>/<workflow_uuid>`.  This is called the `workflow_root` and it is the root directory for all activity in this workflow.

Each `call` has its own subdirectory located at `<workflow_root>/call-<call_name>`.  This is the `<call_dir>`.  For example, having a `stdout` and `stderr` file is common among both backends and they both write a shell script file to the `<call_dir>` as well.  See the descriptions below for details about backend-specific files that are written to these directories.

An example of a workflow output directory for a three-step WDL file might look like this:

```
cromwell-executions/
└── three_step
    └── a59651fc-4d9a-4fed-99ba-f5e2c9d84bb4
        ├── call-cgrep
        │   ├── Users
        │   │   └── jdoe
        │   │       └── projects
        │   │           └── cromwell
        │   │               └── cromwell-executions
        │   │                   └── three_step
        │   │                       └── a59651fc-4d9a-4fed-99ba-f5e2c9d84bb4
        │   │                           └── call-ps
        │   │                               └── stdout
        │   ├── rc
        │   ├── script
        │   ├── stderr
        │   └── stdout
        ├── call-ps
        │   ├── rc
        │   ├── script
        │   ├── stderr
        │   └── stdout
        └── call-wc
            ├── Users
            │   └── jdoe
            │       └── projects
            │           └── cromwell
            │               └── cromwell-executions
            │                   └── three_step
            │                       └── a59651fc-4d9a-4fed-99ba-f5e2c9d84bb4
            │                           └── call-ps
            │                               └── stdout
            ├── rc
            ├── script
            ├── stderr
            └── stdout
```

WDL File

```wdl
task ps {
  command {
    ps
  }
  output {
    File procs = stdout()
  }
}

task cgrep {
  String pattern
  File in_file
  command {
    grep '${pattern}' ${in_file} | wc -l
  }
  output {
    Int count = read_int(stdout())
  }
}

task wc {
  File in_file
  command {
    cat ${in_file} | wc -l
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
```

In the above directory structure, you'll notice that the `call-cgrep` and `call-wc` sub-directories both contain a directory structure to point to the `stdout` file from the invocation of `ps`.  In these cases, that `stdout` file is a localized version of the one within `call-ps/stdout`.  By default both of those `stdout` files would be hard-links but they could also be symbolic links or copies of the file, depending on how Cromwell is configured (see below).  The directory structure is nested so deeply to avoid collisions.  For example, if either of these call invocations referenced two files called `stdout`, they'd collide if they were put into the same directory so the full directory structure is maintained.

Any input files to a call need to be localized into the `<call_dir>`.  There are a few localization strategies that Cromwell will try until one works.  Below is the default order specified in `application.conf` but this order can be overridden:

* `hard-link` - This will create a hard link (not symbolic) link to the file
* `soft-link` - Create a symbolic link to the file.  This strategy is not applicable for tasks which specify a Docker image and will be ignored.
* `copy` - Make a copy the file

Shared filesystem localization is defined in the `config` section of each backend.  The default stanza for the Local, SGE, and associated backends looks like this:

```
filesystems {
 local {
   localization: [
	 "hard-link", "soft-link", "copy"
   ]
 }
}
```

### Google Cloud Storage Filesystem

On the JES backend the GCS (Google Cloud Storage) filesystem is used for the root of the workflow execution.
On the Local, SGE, and associated backends any GCS URI will be downloaded locally.  For the JES backend the `jes_gcs_root` [workflow option](#workflow-options) will take
precedence over the `root` specified at `backend.providers.JES.config.root` in the configuration file. Google Cloud Storage URIs are the only acceptable values for `File` inputs for
workflows using the JES backend.

## Local Backend

The local backend will simply launch a subprocess for each task invocation and wait for it to produce its rc file.

This backend creates three files in the `<call_dir>` (see previous section):

* `script` - A shell script of the job to be run.  This contains the user's command from the `command` section of the WDL code.
* `stdout` - The standard output of the process
* `stderr` - The standard error of the process

The `script` file contains:

```
#!/bin/sh
cd <container_call_root>
<user_command>
echo $? > rc
```

`<container_call_root>` would be equal to `<call_dir>` for non-Docker jobs, or it would be under `/root/<workflow_uuid>/call-<call_name>` if this is running in a Docker container.

When running without docker, the subprocess command that the local backend will launch is:

```
/bin/bash <script>"
```

When running with docker, the subprocess command that the local backend will launch is:

```
docker run --rm -v <cwd>:<docker_cwd> -i <docker_image> /bin/bash < <script>
```

> **NOTE**: If you are using the local backend with Docker and Docker Machine on Mac OS X, by default Cromwell can only
> run from in any path under your home directory.
>
> The `-v` flag will only work if `<cwd>` is within your home directory because VirtualBox with
> Docker Machine only exposes the home directory by default.  Any local path used in `-v` that is not within the user's
> home directory will silently be interpreted as references to paths on the VirtualBox VM.  This can manifest in
> Cromwell as tasks failing for odd reasons (like missing RC file)
>
> See https://docs.docker.com/engine/userguide/dockervolumes/ for more information on volume mounting in Docker.

## Sun GridEngine Backend

The GridEngine and similar backends use programs such as `qsub` to launch a job and will poll the filesystem to determine if a job is completed.

The backend is specified via the actor factory `ConfigBackendLifecycleActorFactory`:

```
backend {
  providers {
    SGE {
      config {
        actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
        # ... other configuration
      }
    }
  }
}
```

This backend makes the same assumption about the filesystem that the local backend does: the Cromwell process and the jobs both have read/write access to the CWD of the job.

The CWD will contain a `script.sh` file which will contain the same contents as the Local backend:

```
#!/bin/sh
cd <container_call_root>
<user_command>
echo $? > rc
```

The job is launched with a configurable command such as:

```bash
qsub \
    -terse \
    -V \
    -b n \
    -N ${job_name} \
    -wd ${cwd} \
    -o ${out} \
    -e ${err} \
    -pe smp ${cpu} \
    ${"-l m_mem_free=" + memory_gb + "gb"} \
    ${"-q " + sge_queue} \
    ${"-P " + sge_project} \
    ${script}
```

The SGE backend gets the job ID from parsing the `submit.stdout` text file.

Since the `script.sh` ends with `echo $? > rc`, the backend will wait for the existence of this file, parse out the return code and determine success or failure and then subsequently post-process.

The command used to submit the job is specified under the configuration key `backend.providers.SGE.config.submit`. It uses the same syntax as a command in WDL, and will be provided the variables:

* `script` - A shell script of the job to be run.  This contains the user's command from the `command` section of the WDL code.
* `cwd` - The path where the script should be run.
* `out` - The path to the stdout.
* `err` - The path to the stderr.
* `job_name` - A unique name for the job.

```
backend {
  providers {
    SGE {
      config {
        # ... other configuration
        submit = """
        qsub \
            -terse \
            -V \
            -b n \
            -N ${job_name} \
            -wd ${cwd} \
            -o ${out} \
            -e ${err} \
            ${script}
        """
      }
    }
  }
}
```

If the backend supports docker, another optional configuration key `backend.providers.<backend>.config.submit-docker` may be specified. When the WDL contains a docker runtime attribute, this command will be provided the two additional variables:

* `docker` - The docker image name.
* `docker_cwd` - The path where `cwd` should be mounted within the docker container.

```
backend {
  providers {
    SGE {
      config {
        # ... other configuration
        submit-docker = """
        qsub \
            -terse \
            -V \
            -b n \
            -N ${job_name} \
            -wd ${cwd} \
            -o ${out} \
            -e ${err} \
            -l docker,docker_images="${docker}"
            -xdv ${cwd}:${docker_cwd}
            ${script}
        """
      }
    }
  }
}
```

If the backend would like to support additional runtime attributes they may be specified in the configuration key `backend.providers.<backend>.config.runtime-attributes`. It uses the same syntax as specifying runtime attributes in a task in WDL.

There are two special runtime attribute configurations, `cpu`, and `memory_<unit>`.

When the runtime attribute configuration `Int cpu` is specified, it is always validated as a positive integer.

When the runtime attribute configuration `Int memory_<unit>` or `Float memory_<unit>` is specified, it is provided to submit by the runtime attribute in WDL `memory`.

For example, if the backend specifies the configuration for `backend.providers.<backend>.config.runtime-attributes` as:

```
backend {
  providers {
    SGE {
      config {
        # ... other configuration
        runtime-attributes = "Float memory_mb"
      }
    }
  }
}
```

And the WDL specifies a task with:

```
task hello_gigabyte {
  command { echo "hello world" }
  runtime { memory: "1 GB" }
}
```

Then for this call, the backend will be provided an additional variable `memory_mb` set to `1000.0`.

Other runtime attributes may be defined by specifying them in under the runtime attributes configuration.

```
backend {
  providers {
    SGE {
      config {
        # ... other configuration
        runtime-attributes = """
        Float memory_mb
        String sge_project
        """
      }
    }
  }
}
```

These variables will then be passed from the WDL into the submit configuration. If one would like to have a default value, just like in WDL, the configuration may specify that the value have a default. The default must match the defined type or an error will be produced.

```
backend {
  providers {
    SGE {
      config {
        # ... other configuration
        runtime-attributes = """
        Float memory_mb = 2.0
        String sge_project = "default"
        """
      }
    }
  }
}
```

Optional values may also be used by appending `?` to the type:

```
backend {
  providers {
    SGE {
      config {
        # ... other configuration
        runtime-attributes = """
        Float? memory_mb
        String? sge_project
        """
      }
    }
  }
}
```

The value will be passed to the submit configuration if provided, and omitted otherwise.

There are also configuration values related to how jobs are rechecked on startup and aborted.

The option is `backend.providers.<backend>.config.run-in-background`. When `true` the backend runs the submit configuration and records the unix process id (PID). To abort the job, the PID is stopped with the unix command `kill`. Upon a cromwell restart, the PID is checked via the unix command `ps` to see if it is still alive, before cromwell goes back to polling for the `rc` file.

When `backend.providers.<backend>.config.run-in-background` is `false`, the default, the backend must specify how read the job identifier from the stdout of the submit, how to kill the job, and how to check if the job is still running during a cromwell restart. These three configuration values are `job-id-regex`, `kill`, and `check-alive`, respectively:

```
backend {
  providers {
    SGE {
      config {
        # ... other configuration
        job-id-regex = "(\\d+)"
        kill = "qdel ${job_id}"
        check-alive = "qstat -j ${job_id}"
        """
      }
    }
  }
}
```

The `job-id-regex` should contain one capture group while matching against the whole line or stdout file. The `check-alive` should return zero if the job is still alive.

## HtCondor Backend

Allows to execute jobs using HTCondor which is a specialized workload management system for compute-intensive jobs created by the Center for High Throughput Computing in the Department of Computer Sciences at the University of Wisconsin-Madison (UW-Madison).

This backend creates six files in the `<call_dir>` (see previous section):

* `script` - A shell script of the job to be run.  This contains the user's command from the `command` section of the WDL code.
* `stdout` - The standard output of the process
* `stderr` - The standard error of the process
* `submitfile` - A submit file that HtCondor understands in order to submit a job
* `submitfile.stdout` - The standard output of the submit file
* `submitfile.stderr` - The standard error of the submit file

The `script` file contains:

```
cd <container_call_root>
<user_command>
echo $? > rc
```

The `submitfile` file contains:

```
executable=cromwell-executions/test/e950e07d-4132-4fe0-8d86-ab6925dd94ad/call-merge_files/script
output=cromwell-executions/test/e950e07d-4132-4fe0-8d86-ab6925dd94ad/call-merge_files/stdout
error=cromwell-executions/test/e950e07d-4132-4fe0-8d86-ab6925dd94ad/call-merge_files/stderr
log=cromwell-executions/test/e950e07d-4132-4fe0-8d86-ab6925dd94ad/call-merge_files/merge_files.log
queue

```

### Caching configuration
This implementation also add basic caching support. It relies in a cache provider to store successful job results.
By default a MongoDB based cache implementation is provided but there is the option of implementing a new provider based on CacheActorFactory and CacheActor interfaces.

From application.conf file:
```
cache {
  provider = "cromwell.backend.impl.htcondor.caching.provider.mongodb.MongoCacheActorFactory"
  enabled = true
  forceRewrite = false
  db {
    host = "127.0.0.1"
    port = 27017
    name = "htcondor"
    collection = "cache"
  }
}

```

* provider: it defines the provider to use based on CacheActorFactory and CacheActor interfaces.
* enabled: enables or disables cache.
* forceRewrite: it allows to invalidate the cache entry and store result again.
* db section: configuration related to MongoDB provider. It may not exist for other implementations.

### Docker
This backend supports the following optional runtime attributes / workflow options for working with Docker:
* docker: Docker image to use such as "Ubuntu".
* dockerWorkingDir: defines the working directory in the container.
* dockerOutputDir: defiles the output directory in the container when there is the need to define a volume for outputs within the container. By default if this attribute is not set, dockerOutputDir will be the job working directory.

Inputs:
HtCondor backend analyzes all inputs and do a distinct of the folders in order to mount input folders into the container.

Outputs:
It will use dockerOutputDir runtime attribute / workflow option to resolve the folder in which the execution results will placed. If there is no dockerOutputDir defined it will use the current working directory.

### CPU, Memory and Disk
This backend supports CPU, memory and disk size configuration through the use of the following runtime attributes / workflow options:
* cpu: defines the amount of CPU to use. Default value: 1. Type: Integer. Ex: 4.
* memory: defines the amount of memory to use. Default value: "512 MB". Type: String. Ex: "4 GB" or "4096 MB"
* disk: defines the amount of disk to use. Default value: "1024 MB". Type: String. Ex: "1 GB" or "1024 MB"

It they are not set, HtCondor backend will use default values.

### Native Specifications
The use of runtime attribute 'nativeSpecs' allows to the user to attach custom HtCondor configuration to tasks.
An example of this is when there is a need to work with 'requirements' or 'rank' configuration.

```
"runtimeAttributes": {
    "nativeSpecs": ["requirements = Arch == \"INTEL\"", "rank = Memory >= 64"]
}
```

nativeSpecs attribute needs to be specified as an array of strings to work.

## Spark Backend

This backend adds support for execution of spark jobs in a workflow using the existing wdl format. 

It supports the following Spark deploy modes:

*  Client deploy mode using the spark standalone cluster manager
*  Cluster deploy mode using the spark standalone cluster manager
*  Client deploy mode using Yarn resource manager
*  Cluster deploy mode using Yarn resource manager

### Configuring Spark Project

When using Spark backend uncomment the following Spark configuration in the application.conf file

```conf
Spark {
       actor-factory = "cromwell.backend.impl.spark.SparkBackendFactory"
       config {
         # Root directory where Cromwell writes job results.  This directory must be
         # visible and writeable by the Cromwell process as well as the jobs that Cromwell
         # launches.
         root: "cromwell-executions"

         filesystems {
           local {
             localization: [
               "hard-link", "soft-link", "copy"
             ]
           }
         }
		master: "local"
		deployMode: "client"
        }

      }
```
and add backend provider as Spark. 

```
backend {
  default = "Spark"
  providers {
  ....
```

### Configuring Spark Master and Deploy Mode

Default configuration is as follows:

```conf
Spark {
		......
		master: "local"
		deployMode: "client"

      }
```

However to use Spark in standalone cluster mode change `master: spark://hostname:6066` and `deployMode: cluster` similarly, for yarn change `master: yarn` and `deployMode: cluster` or `deployMode: client` to run in cluster or client mode respectively. 

### Spark runtime attributes

Supported runtime attributes for a Spark Job is as follows:
	
* executorCores (default value is 1)
* executorMemory (default value is "1 GB")
* appMainClass ( Spark app/job entry point)
* numberOfExecutors ( Specific to cluster deploy mode)

Sample usage :

```wdl
task sparkjob_with_yarn_cluster {
        .....
        
        runtime {
                appMainClass: "${entry_point}"
                executorMemory: "4G"
                executorCores: "2"
        }
        
        .....
	}
```

### Spark Environment

The Spark backend assumes Spark is already installed, and it constructs the spark submit command with the `SPARK_HOME` environment variable if set. Otherwise backend creates command `spark-submit` without a fully qualified path to `spark-submit`.

Supported File Systems as follows: 

* Local File System
* Network File System
* Distributed file system

### Sample Wdl
Next, create a Wdl, and it's json input like so:

```wdl
task sparkjob_with_yarn_cluster {
        File input_jar
        String input_1
        String output_base
        String entry_point
        Int cores
        String memory

        command {
                ${input_jar} ${input_1} ${output_base}
        }

        runtime {
                appMainClass: "${entry_point}"
                executorMemory: "${memory}"
                executorCores: "${cores}"
        }
        output {
                File out = "${output_base}"
          }
	}
```

and its accompanying json input as:

```json
{
	"sparkWithYarnCluster.sparkjob_with_yarn_cluster.memory": "4G",
	"sparkWithYarnCluster.sparkjob_with_yarn_cluster.output_base":"/mnt/lustre/hadoop/home/yarn_cluster_output",
	"sparkWithYarnCluster.sparkjob_with_yarn_cluster.entry_point": "com.org.spark.poc.nfs.SparkVowelLine",
	"sparkWithYarnCluster.sparkjob_with_yarn_cluster.cores": "12",
	"sparkWithYarnCluster.sparkjob_with_yarn_cluster.input_1": "/mnt/lustre/hadoop/home/inputfiles/sample.txt",
	"sparkWithYarnCluster.sparkjob_with_yarn_cluster.input_jar": "/mnt/lustre/hadoop/home/inputjars/spark_hdfs.jar"
}
```

## Google JES Backend

Google JES (Job Execution Service) is a Docker-as-a-service from Google.

### Configuring Google Project

You'll need the following things to get started:

* A Google Project (Manage/create projects [here](https://console.developers.google.com/project))
* A Google Cloud Storage bucket (View/create buckets in your project [here](https://console.cloud.google.com/storage/browser))

On your Google project, open up the [API Manager](https://console.developers.google.com/apis/library) and enable the following APIs:

* Google Compute Engine
* Google Cloud Storage
* Genomics API

If your project is `my-project` your bucket is `gs://my-bucket/`, then update your [Cromwell configuration file](#configuring-cromwell) as follows:

```hocon
backend {
  default = "JES"
  providers {
    JES {
      actor-factory = "cromwell.backend.impl.jes.JesBackendLifecycleActorFactory"
      config {
        project = "my-project"
        root = "gs://my-bucket"
        genomics-api-queries-per-100-seconds = 1000
        .
        .
        .
      }
    }
  ]
}
```

If your project has API quotas other than the defaults set the `genomics-api-queries-per-100-seconds` value to be the lesser of the `Queries per 100 seconds per user` and `Queries per 100 seconds` quotas. This value will be used to help tune Cromwell's rate of interaction with JES.

### Configuring Authentication

The `google` stanza in the Cromwell configuration file defines how to authenticate to Google.  There are four different
authentication schemes that might be used:

* `application_default` - (default, recommended) Use [application default](https://developers.google.com/identity/protocols/application-default-credentials) credentials.
* `service_account` - Use a specific service account and key file (in PEM format) to authenticate.
* `user_account` - Authenticate as a user.
* `refresh_token` - Authenticate using a refresh token supplied in the workflow options.

The `auths` block in the `google` stanza defines the authorization schemes within a Cromwell deployment:

```hocon
google {
  application-name = "cromwell"
  auths = [
    {
      name = "application-default"
      scheme = "application_default"
    },
    {
      name = "user-via-refresh"
      scheme = "refresh_token"
      client-id = "secret_id"
      client-secret = "secret_secret"
    },
    {
      name = "service-account"
      scheme = "service_account"
      service-account-id = "my-service-account"
      pem-file = "/path/to/file.pem"
    }
  ]
}
```

These authorization schemes can be referenced by name within other portions of the configuration file.  For example, both
the `genomics` and `filesystems.gcs` sections within a JES configuration block must reference an auth defined in this block.
The auth for the `genomics` section governs the interactions with JES itself, while `filesystems.gcs` governs the localization
of data into and out of GCE VMs.

#### Application Default Credentials

By default, application default credentials will be used.  There is no configuration required for application default
credentials, only `name` and `scheme` are required.

To authenticate, run the following commands from your command line (requires [gcloud](https://cloud.google.com/sdk/gcloud/)):

```
$ gcloud auth login
$ gcloud config set project my-project
```

#### Service Account

First create a new service account through the [API Credentials](https://console.developers.google.com/apis/credentials) page.  Go to **Create credentials -> Service account key**.  Then in the **Service account** dropdown select **New service account**.  Fill in a name (e.g. `my-account`), and select key type of JSON.

Creating the account will cause the JSON file to be downloaded.  The structure of this file is roughly like this (account name is `my-account`):

```
{
  "type": "service_account",
  "project_id": "my-project",
  "private_key_id": "OMITTED",
  "private_key": "-----BEGIN PRIVATE KEY-----\nBASE64 ENCODED KEY WITH \n TO REPRESENT NEWLINES\n-----END PRIVATE KEY-----\n",
  "client_email": "my-account@my-project.iam.gserviceaccount.com",
  "client_id": "22377410244549202395",
  "auth_uri": "https://accounts.google.com/o/oauth2/auth",
  "token_uri": "https://accounts.google.com/o/oauth2/token",
  "auth_provider_x509_cert_url": "https://www.googleapis.com/oauth2/v1/certs",
  "client_x509_cert_url": "https://www.googleapis.com/robot/v1/metadata/x509/my-account%40my-project.iam.gserviceaccount.com"
}
```

Most importantly, the value of the `client_email` field should go into the `service-account-id` field in the configuration (see below).  The
`private_key` portion needs to be pulled into its own file (e.g. `my-key.pem`).  The `\n`s in the string need to be converted to newline characters.

While technically not part of Service Account authorization mode, one can also override the default service account that the compute VM is started with via the configuration option `JES.config.genomics.compute-service-account` or through the workflow options parameter `google_compute_service_account`.  It's important that this service account, and the service account specified in `JES.config.genomics.auth` can both read/write the location specified by `JES.config.root`

#### Refresh Token

A **refresh_token** field must be specified in the [workflow options](#workflow-options) when submitting the job.  Omitting this field will cause the workflow to fail.

The refresh token is passed to JES along with the `client-id` and `client-secret` pair specified in the corresponding entry in `auths`.


### Docker

It is possible to reference private docker images in DockerHub to be run on JES.
However, in order for the image to be pulled, the docker credentials with access to this image must be provided in the configuration file.


```
backend {
  default = "JES"
  providers {
    JES {
      actor-factory = "cromwell.backend.impl.local.LocalBackendLifecycleActorFactory"
      config {
        dockerhub {
          account = "mydockeraccount@mail.com"
          token = "mydockertoken"
        }
      }
    }
  }
}
```

It is now possible to reference an image only this account has access to:

```
task mytask {
  command {
    ...
  }
  runtime {
    docker: "private_repo/image"
    memory: "8 GB"
    cpu: "1"
  }
  ...
}
```

Note that if the docker image to be used is public there is no need to add this configuration.

### Monitoring

In order to monitor metrics (CPU, Memory, Disk usage...) about the VM during Call Runtime, a workflow option can be used to specify the path to a script that will run in the background and write its output to a log file.

```
{
  "monitoring_script": "gs://cromwell/monitoring/script.sh"
}
```

The output of this script will be written to a `monitoring.log` file that will be available in the call gcs bucket when the call completes.  This feature is meant to run a script in the background during long-running processes.  It's possible that if the task is very short that the log file does not flush before de-localization happens and you will end up with a zero byte file.

# Runtime Attributes

Runtime attributes are used to customize tasks. Within a task one can specify runtime attributes to customize the environment for the call.

For example:

```
task jes_task {
  command {
    echo "Hello JES!"
  }
  runtime {
    docker: "ubuntu:latest"
    memory: "4G"
    cpu: "3"
    zones: "us-central1-c us-central1-b"
    disks: "/mnt/mnt1 3 SSD, /mnt/mnt2 500 HDD"
  }
}
workflow jes_workflow {
  call jes_task
}
```

This table lists the currently available runtime attributes for cromwell:

| Runtime Attribute    | LOCAL |  JES  |  SGE  |
| -------------------- |:-----:|:-----:|:-----:|
| continueOnReturnCode |   x   |   x   |   x   |
| cpu                  |       |   x   |   x   |
| disks                |       |   x   |       |
| zones                |       |   x   |       |
| docker               |   x   |   x   |   x   |
| failOnStderr         |   x   |   x   |   x   |
| memory               |       |   x   |   x   |
| preemptible          |       |   x   |       |
| bootDiskSizeGb       |       |   x   |       |

Runtime attribute values are interpreted as expressions.  This means that it is possible to express the value of a runtime attribute as a function of one of the task's inputs.  For example:

```
task runtime_test {
  String ubuntu_tag
  Int memory_gb

  command {
    ./my_binary
  }

  runtime {
    docker: "ubuntu:" + ubuntu_tag
    memory: memory_gb + "GB"
  }
}
```

SGE and similar backends may define other configurable runtime attributes beyond the five listed. See [Sun GridEngine](#sun-gridengine-backend) for more information.

## Specifying Default Values

Default values for runtime attributes can be specified via [workflow options](#workflow-options).  For example, consider this WDL file:

```wdl
task first {
  command { ... }
}

task second {
  command {...}
  runtime {
    docker: "my_docker_image"
  }
}

workflow w {
  call first
  call second
}
```

And this set of workflow options:

```json
{
  "default_runtime_attributes": {
    "docker": "ubuntu:latest",
    "zones": "us-central1-c us-central1-b"
  }
}
```

Then these values for `docker` and `zones` will be used for any task that does not explicitly override them in the WDL file. So the effective runtime for `task first` is:
```
{
    "docker": "ubuntu:latest",
    "zones": "us-central1-c us-central1-b"
  }
```
And the effective runtime for `task second` is:
```
{
    "docker": "my_docker_image",
    "zones": "us-central1-c us-central1-b"
  }
```
Note how for task second, the WDL value for `docker` is used instead of the default provided in the workflow options.

## continueOnReturnCode

When each task finishes it returns a code. Normally, a non-zero return code indicates a failure. However you can override this behavior by specifying the `continueOnReturnCode` attribute.

When set to false, any non-zero return code will be considered a failure. When set to true, all return codes will be considered successful.

```
runtime {
  continueOnReturnCode: true
}
```

When set to an integer, or an array of integers, only those integers will be considered as successful return codes.

```
runtime {
  continueOnReturnCode: 1
}
```

```
runtime {
  continueOnReturnCode: [0, 1]
}
```

Defaults to "0".

## cpu

Passed to JES: "The minimum number of cores to use."

Passed to SGE, etc.: Configurable, but usually a reservation and/or limit of number of cores.

```
runtime {
  cpu: 2
}
```

Defaults to "1".

## disks

Passed to JES: "Disks to attach."

The disks are specified as a comma separated list of disks. Each disk is further separated as a space separated triplet of:

1. Mount point (absolute path), or `local-disk` to reference the mount point where JES will localize files and the task's current working directory will be
2. Disk size in GB (ignored for disk type LOCAL)
3. Disk type.  One of: "LOCAL", "SSD", or "HDD" ([documentation](https://cloud.google.com/compute/docs/disks/#overview))

All tasks launched on JES *must* have a `local-disk`.  If one is not specified in the runtime section of the task, then a default of `local-disk 10 SSD` will be used.  The `local-disk` will be mounted to `/cromwell_root`.

The Disk type must be one of "LOCAL", "SSD", or "HDD". When set to "LOCAL", the size of the drive is automatically provisioned by Google so any size specified in WDL will be ignored. All disks are set to auto-delete after the job completes.

**Example 1: Changing the Localization Disk**

```
runtime {
  disks: "local-disk 100 SSD"
}
```

**Example 2: Mounting an Additional Two Disks**

```
runtime {
  disks: "/mnt/my_mnt 3 SSD, /mnt/my_mnt2 500 HDD"
}
```

### Boot Disk
In addition to working disks, JES allows specification of a boot disk size. This is the disk where the docker image itself is booted, **not the working directory of your task on the VM**.
Its primary purpose is to ensure that larger docker images can fit on the boot disk.
```
runtime {
  # Yikes, we have a big OS in this docker image! Allow 50GB to hold it:
  bootDiskSizeGb: 50
}
```

Since no `local-disk` entry is specified, Cromwell will automatically add `local-disk 10 SSD` to this list.

## zones

The ordered list of zone preference (see [Region and Zones](https://cloud.google.com/compute/docs/zones) documentation for specifics)

The zones are specified as a space separated list, with no commas.

```
runtime {
  zones: "us-central1-c us-central1-b"
}
```

Defaults to "us-central1-b"

## docker

When specified, cromwell will run your task within the specified Docker image.

```
runtime {
  docker: "ubuntu:latest"
}
```

This attribute is mandatory when submitting tasks to JES. When running on other backends, they default to not running the process within Docker.

## failOnStderr

Some programs write to the standard error stream when there is an error, but still return a zero exit code. Set `failOnStderr` to true for these tasks, and it will be considered a failure if anything is written to the standard error stream.

```
runtime {
  failOnStderr: true
}
```

Defaults to "false".

## memory

Passed to JES: "The minimum amount of RAM to use."

Passed to SGE, etc.: Configurable, but usually a reservation and/or limit of memory.

The memory size is specified as an amount and units of memory, for example "4 G".

```
runtime {
  memory: "4G"
}
```

Defaults to "2G".

## preemptible

Passed to JES: "If applicable, preemptible machines may be used for the run."

Take an Int as a value that indicates the maximum number of times Cromwell should request a preemptible machine for this task before defaulting back to a non-preemptible one.
eg. With a value of 1, Cromwell will request a preemptible VM, if the VM is preempted, the task will be retried with a non-preemptible VM.
Note: If specified, this attribute overrides [workflow options](#workflow-options).

```
runtime {
  preemptible: 1
}
```

Defaults to "false".

# Logging

Cromwell accepts two Java Properties or Environment Variables for controlling logging:

* `LOG_MODE` - Accepts either `pretty` or `standard` (default `pretty`).  In `standard` mode, logs will be written without ANSI escape code coloring, with a layout more appropriate for server logs, versus `pretty` that is easier to read for a single workflow run.
* `LOG_LEVEL` - Level at which to log (default `info`).

Additionally, a directory may be set for writing per workflow logs. By default, the per workflow logs will be erased once the workflow completes.

```hocon
// In application.conf or specified via system properties
workflow-options {
    workflow-log-dir: "cromwell-workflow-logs"
    workflow-log-temporary: true
}
```

The usual case of generating the temporary per workflow logs is to copy them to a remote directory, while deleting the local copy to preserve local disk space. To specify the remote directory to copy the logs to use the separate [workflow option](#workflow-options) `final_workflow_log_dir`.

# Workflow Options

When running a workflow from the [command line](#run) or [REST API](#post-apiworkflowsversion), one may specify a JSON file that toggles various options for running the workflow.  From the command line, the workflow options is passed in as the third positional parameter to the 'run' subcommand.  From the REST API, it's an optional part in the multi-part POST request.  See the respective sections for more details.

Example workflow options file:

```json
{
  "jes_gcs_root": "gs://my-bucket/workflows",
  "google_project": "my_google_project",
  "refresh_token": "1/Fjf8gfJr5fdfNf9dk26fdn23FDm4x"
}
```

Valid keys and their meanings:

* Global *(use with any backend)*
    * **write_to_cache** - Accepts values `true` or `false`.  If `false`, the completed calls from this workflow will not be added to the cache.  See the [Call Caching](#call-caching) section for more details.
    * **read_from_cache** - Accepts values `true` or `false`.  If `false`, Cromwell will not search the cache when invoking a call (i.e. every call will be executed unconditionally).  See the [Call Caching](#call-caching) section for more details.
    * **final_workflow_log_dir** - Specifies a path where per-workflow logs will be written.  If this is not specified, per-workflow logs will not be copied out of the Cromwell workflow log temporary directory/path before they are deleted.
    * **final_workflow_outputs_dir** - Specifies a path where final workflow outputs will be written.  If this is not specified, workflow outputs will not be copied out of the Cromwell workflow execution directory/path.
    * **final_call_logs_dir** - Specifies a path where final call logs will be written.  If this is not specified, call logs will not be copied out of the Cromwell workflow execution directory/path.
    * **default_runtime_attributes** - A JSON object where the keys are [runtime attributes](#runtime-attributes) and the values are defaults that will be used through the workflow invocation.  Individual tasks can choose to override these values.  See the [runtime attributes](#specifying-default-values) section for more information.
    * **continueOnReturnCode** - Can accept a boolean value or a comma separated list of integers in a string.  Defaults to false.  If false, then only return code of 0 will be acceptable for a task invocation.  If true, then any return code is valid.  If the value is a list of comma-separated integers in a string, this is interpreted as the acceptable return codes for this task.
    * **workflow_failure_mode** - What happens after a task fails. Choose from:
        * **ContinueWhilePossible** - continues to start and process calls in the workflow, as long as they did not depend on the failing call
        * **NoNewCalls** - no *new* calls are started but existing calls are allowed to finish
        * The default is `NoNewCalls` but this can be changed using the `workflow-options.workflow-failure-mode` configuration option.
    * **backend** - Override the default backend specified in the Cromwell configuration for this workflow only.
* JES Backend Only
    * **jes_gcs_root** - (JES backend only) Specifies where outputs of the workflow will be written.  Expects this to be a GCS URL (e.g. `gs://my-bucket/workflows`).  If this is not set, this defaults to the value within `backend.jes.config.root` in the [configuration](#configuring-cromwell).
    * **google_compute_service_account** - (JES backend only) Specifies an alternate service account to use on the compute instance (e.g. my-new-svcacct@my-google-project.iam.gserviceaccount.com).  If this is not set, this defaults to the value within `backend.jes.config.genomics.compute-service-account` in the [configuration](#configuring-cromwell) if specified or `default` otherwise.
    * **google_project** - (JES backend only) Specifies which google project to execute this workflow.
    * **refresh_token** - (JES backend only) Only used if `localizeWithRefreshToken` is specified in the [configuration file](#configuring-cromwell).
    * **auth_bucket** - (JES backend only) defaults to the the value in **jes_gcs_root**.  This should represent a GCS URL that only Cromwell can write to.  The Cromwell account is determined by the `google.authScheme` (and the corresponding `google.userAuth` and `google.serviceAuth`)
    * **monitoring_script** - (JES backend only) Specifies a GCS URL to a script that will be invoked prior to the WDL command being run.  For example, if the value for monitoring_script is "gs://bucket/script.sh", it will be invoked as `./script.sh > monitoring.log &`.  The value `monitoring.log` file will be automatically de-localized.

# Call Caching

Call Caching allows Cromwell to detect when a job has been run in the past so it doesn't have to re-compute results.  Cromwell searches the cache of previously run jobs for a one that has the exact same command and exact same inputs.  If a previously run job is found in the cache, Cromwell will **copy the results** of the previous job instead of re-running it.

Cromwell's call cache is maintained in its database.  For best mileage with call caching, configure Cromwell to [point to a MySQL database](#database) instead of the default in-memory database.  This way any invocation of Cromwell (either with `run` or `server` subcommands) will be able to utilize results from all calls that are in that database.

**Call Caching is disabled by default.**  Once enabled, Cromwell will search the call cache for every `call` statement invocation, assuming `read_from_cache` is enabled (see below):

* If there was no cache hit, the `call` will be executed as normal.  Once finished it will add itself to the cache, assuming `read_from_cache` is enabled (see below)
* If there was a cache hit, outputs are **copied** from the cached job to the new job's output directory

> **Note:** If call caching is enabled, be careful not to change the contents of the output directory for any previously run job.  Doing so might cause cache hits in Cromwell to copy over modified data and Cromwell currently does not check that the contents of the output directory changed.

## Configuring Call Caching
To enable Call Caching, add the following to your Cromwell [configuration](#configuring-cromwell):

```
call-caching {
  enabled = true
  invalidate-bad-cache-results = true
}
```

When `call-caching.enabled=true` (default: `false`), Cromwell will be able to to copy results from previously run jobs (when appropriate).
When `invalidate-bad-cache-results=true` (default: `true`), Cromwell will invalidate any cache results which fail to copy during a cache-hit. This is usually desired but might be unwanted if a cache might fail to copy for external reasons, such as a difference in user authentication.

## Call Caching Workflow Options
Cromwell also accepts two [workflow option](#workflow-options) related to call caching:

* If call caching is enabled, but one wishes to run a workflow but not add any of the calls into the call cache when they finish, the `write_to_cache` option can be set to `false`.  This value defaults to `true`.
* If call caching is enabled, but you don't want to check the cache for any `call` invocations, set the option `read_from_cache` to `false`.  This value also defaults to `true`

> **Note:** If call caching is disabled, the to workflow options `read_from_cache` and `write_to_cache` will be ignored and the options will be treated as though they were 'false'.

## Local Filesystem Options
When running a job on the Config (Shared Filesystem) backend, Cromwell provides some additional options in the backend's config section:

```
      config {
        ...
        filesystems {
          ...
          local {
            ...
            caching {
              # When copying a cached result, what type of file duplication should occur. Attempted in the order listed below:
              duplication-strategy: [
                "hard-link", "soft-link", "copy"
              ]

              # Possible values: file, path
              # "file" will compute an md5 hash of the file content.
              # "path" will compute an md5 hash of the file path. This strategy will only be effective if the duplication-strategy (above) is set to "soft-link",
              # in order to allow for the original file path to be hashed.
              # Default: file
              hashing-strategy: "file"

              # When true, will check if a sibling file with the same name and the .md5 extension exists, and if it does, use the content of this file as a hash.
              # If false or the md5 does not exist, will proceed with the above-defined hashing strategy.
              # Default: false
              check-sibling-md5: false
            }
          }
        }
      }
```
# Imports

Import statements inside of a WDL file are supported by Cromwell when running in Server mode as well as Single Workflow Runner Mode.

In Single Workflow Runner Mode, you pass in a zip file which includes the WDL files referenced by the import statements. Cromwell requires the zip file to be passed in as a command line argument, as explained by the section [run](#run).

For example, given a workflow `wf.wdl` and an imports directory `WdlImports.zip`, a sample command would be:
```
java -jar cromwell.jar wf.wdl wf.inputs - - WdlImports.zip
```

In Server Mode, you pass in a zip file using the parameter `wdlDependencies` via the [POST /api/workflows/:version](#post-apiworkflowsversion) endpoint.


# Sub Workflows

WDL allows the execution of an entire workflow as a step in a larger workflow (see WDL SPEC for more details), which is what will be referred to as a sub workflow going forward.
Cromwell supports execution of such workflows. Note that sub workflows can themselves contain sub workflows, etc... There is no limitation as to how deeply workflows can be nested.

## Execution

Sub workflows are executed exactly as a task would be.
*This means that if another call depends on an output of a sub workflow, this call will run when the whole sub workflow completes (successfully).*
For example, in the following case :

`main.wdl`
```
import "sub_wdl.wdl" as sub

workflow main_workflow {

    call sub.hello_and_goodbye { input: hello_and_goodbye_input = "sub world" }
    
    # call myTask { input: hello_and_goodbye.hello_output }
    
    output {
        String main_output = hello_and_goodbye.hello_output
    }
}
```

`sub_wdl.wdl`
```
task hello {
  String addressee
  command {
    echo "Hello ${addressee}!"
  }
  runtime {
      docker: "ubuntu:latest"
  }
  output {
    String salutation = read_string(stdout())
  }
}

task goodbye {
  String addressee
  command {
    echo "Goodbye ${addressee}!"
  }
  runtime {
      docker: "ubuntu:latest"
  }
  output {
    String salutation = read_string(stdout())
  }
}

workflow hello_and_goodbye {
  String hello_and_goodbye_input
  
  call hello {input: addressee = hello_and_goodbye_input }
  call goodbye {input: addressee = hello_and_goodbye_input }
  
  output {
    String hello_output = hello.salutation
    String goodbye_output = goodbye.salutation
  }
}
```

`myTask` will start only when hello_and_goodbye completes (which means all of its calls are done), even though `myTask` only needs the output of hello in the hello_and_goodbye sub workflow. 
If hello_and_goodbye fails, then `myTask` won't be executed.
Only workflow outputs are visible outside a workflow, which means that references to outputs produced by a sub workflow will only be valid if those outputs are exposed in the workflow output section.

Sub workflows are executed in the context of a main workflow, which means that operations that are normally executed once per workflow (set up, clean up, outputs copying, log copying, etc...)
will NOT be re-executed for each sub workflow. For instance if a resource is created during workflow initialization, sub workflows will need to share this same resource.
Workflow outputs will be copied for the main root workflow but not for intermediate sub workflows.

Restarts, aborts, and call-caching work exactly as they would with tasks. 
All tasks run by a sub workflow are eligible for call caching under the same rules as any other task.
However, workflows themselves are not cached as such. Which means that running the exact same workflow twice with call caching on will trigger each task to cache individually,
but not the workflow itself.

The root path for sub workflow execution files (scripts, output files, logs) will be under the parent workflow call directory.
For example, the execution directory for the above main workflow would look like the following:

```
cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/ <- main workflow id
└── call-hello_and_goodbye <- call directory for call hello_and_goodbye in the main workflow
    └── hello_and_goodbye <- name of the sub workflow 
        └── a6365f91-c807-465a-9186-a5d3da98fe11 <- sub workflow id
            ├── call-goodbye
            │   └── execution
            │       ├── rc
            │       ├── script
            │       ├── script.background
            │       ├── script.submit
            │       ├── stderr
            │       ├── stderr.background
            │       ├── stdout
            │       └── stdout.background
            └── call-hello
                └── execution
                    ├── rc
                    ├── script
                    ├── script.background
                    ├── script.submit
                    ├── stderr
                    ├── stderr.background
                    ├── stdout
                    └── stdout.background

```

## Metadata
Each sub workflow will have its own workflow ID. This ID will appear in the metadata of the parent workflow, in the call section corresponding to the sub workflow, under the "subWorkflowId" attribute.
For example, querying the `main_workflow` metadata above (minus the `myTask` call) , could result in something like this:

`GET /api/workflows/v2/1d919bd4-d046-43b0-9918-9964509689dd/metadata`

```
{
  "workflowName": "main_workflow",
  "submittedFiles": {
    "inputs": "{}",
    "workflow": "import \"sub_wdl.wdl\" as sub\n\nworkflow main_workflow {\n\n    call sub.hello_and_goodbye { input: hello_and_goodbye_input = \"sub world\" }\n    \n    # call myTask { input: hello_and_goodbye.hello_output }\n    \n    output {\n        String main_output = hello_and_goodbye.hello_output\n    }\n}",
    "options": "{\n\n}"
  },
  "calls": {
    "main_workflow.hello_and_goodbye": [
      {
        "executionStatus": "Done",
        "shardIndex": -1,
        "outputs": {
          "goodbye_output": "Goodbye sub world!",
          "hello_output": "Hello sub world!"
        },
        "inputs": {
          "hello_and_goodbye_input": "sub world"
        },
        "end": "2016-11-17T14:13:41.117-05:00",
        "attempt": 1,
        "start": "2016-11-17T14:13:39.236-05:00",
        "subWorkflowId": "a6365f91-c807-465a-9186-a5d3da98fe11"
      }
    ]
  },
  "outputs": {
    "main_output": "Hello sub world!"
  },
  "workflowRoot": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd",
  "id": "1d919bd4-d046-43b0-9918-9964509689dd",
  "inputs": {},
  "submission": "2016-11-17T14:13:39.104-05:00",
  "status": "Succeeded",
  "end": "2016-11-17T14:13:41.120-05:00",
  "start": "2016-11-17T14:13:39.204-05:00"
}
```

The sub workflow ID can be queried separately:

`GET /api/workflows/v2/a6365f91-c807-465a-9186-a5d3da98fe11/metadata`

```
{
  "workflowName": "hello_and_goodbye",
  "calls": {
    "sub.hello_and_goodbye.hello": [
      {
        "executionStatus": "Done",
        "stdout": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-hello/execution/stdout",
        "shardIndex": -1,
        "outputs": {
          "salutation": "Hello sub world!"
        },
        "runtimeAttributes": {
          "docker": "ubuntu:latest",
          "failOnStderr": false,
          "continueOnReturnCode": "0"
        },
        "cache": {
          "allowResultReuse": true
        },
        "Effective call caching mode": "CallCachingOff",
        "inputs": {
          "addressee": "sub world"
        },
        "returnCode": 0,
        "jobId": "49830",
        "backend": "Local",
        "end": "2016-11-17T14:13:40.712-05:00",
        "stderr": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-hello/execution/stderr",
        "callRoot": "/cromwell/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-hello",
        "attempt": 1,
        "executionEvents": [
          {
            "startTime": "2016-11-17T14:13:39.240-05:00",
            "description": "Pending",
            "endTime": "2016-11-17T14:13:39.240-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:39.240-05:00",
            "description": "RequestingExecutionToken",
            "endTime": "2016-11-17T14:13:39.240-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:39.240-05:00",
            "description": "PreparingJob",
            "endTime": "2016-11-17T14:13:39.243-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:39.243-05:00",
            "description": "RunningJob",
            "endTime": "2016-11-17T14:13:40.704-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:40.704-05:00",
            "description": "UpdatingJobStore",
            "endTime": "2016-11-17T14:13:40.712-05:00"
          }
        ],
        "start": "2016-11-17T14:13:39.239-05:00"
      }
    ],
    "sub.hello_and_goodbye.goodbye": [
      {
        "executionStatus": "Done",
        "stdout": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-goodbye/execution/stdout",
        "shardIndex": -1,
        "outputs": {
          "salutation": "Goodbye sub world!"
        },
        "runtimeAttributes": {
          "docker": "ubuntu:latest",
          "failOnStderr": false,
          "continueOnReturnCode": "0"
        },
        "cache": {
          "allowResultReuse": true
        },
        "Effective call caching mode": "CallCachingOff",
        "inputs": {
          "addressee": "sub world"
        },
        "returnCode": 0,
        "jobId": "49831",
        "backend": "Local",
        "end": "2016-11-17T14:13:41.115-05:00",
        "stderr": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-goodbye/execution/stderr",
        "callRoot": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-goodbye",
        "attempt": 1,
        "executionEvents": [
          {
            "startTime": "2016-11-17T14:13:39.240-05:00",
            "description": "Pending",
            "endTime": "2016-11-17T14:13:39.240-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:39.240-05:00",
            "description": "RequestingExecutionToken",
            "endTime": "2016-11-17T14:13:39.240-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:39.240-05:00",
            "description": "PreparingJob",
            "endTime": "2016-11-17T14:13:39.243-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:39.243-05:00",
            "description": "RunningJob",
            "endTime": "2016-11-17T14:13:41.112-05:00"
          },
          {
            "startTime": "2016-11-17T14:13:41.112-05:00",
            "description": "UpdatingJobStore",
            "endTime": "2016-11-17T14:13:41.115-05:00"
          }
        ],
        "start": "2016-11-17T14:13:39.239-05:00"
      }
    ]
  },
  "outputs": {
    "goodbye_output": "Goodbye sub world!",
    "hello_output": "Hello sub world!"
  },
  "workflowRoot": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11",
  "id": "a6365f91-c807-465a-9186-a5d3da98fe11",
  "inputs": {
    "hello_and_goodbye_input": "sub world"
  },
  "status": "Succeeded",
  "parentWorkflowId": "1d919bd4-d046-43b0-9918-9964509689dd",
  "end": "2016-11-17T14:13:41.116-05:00",
  "start": "2016-11-17T14:13:39.236-05:00"
}
```

It's also possible to set the URL query parameter `expandSubWorkflows` to `true` to automatically include sub workflows metadata (`false` by default).

`GET api/workflows/v2/1d919bd4-d046-43b0-9918-9964509689dd/metadata?expandSubWorkflows=true`

```
{
  "workflowName": "main_workflow",
  "submittedFiles": {
    "inputs": "{}",
    "workflow": "import \"sub_wdl.wdl\" as sub\n\nworkflow main_workflow {\n\n    call sub.hello_and_goodbye { input: hello_and_goodbye_input = \"sub world\" }\n    \n    # call myTask { input: hello_and_goodbye.hello_output }\n    \n    output {\n        String main_output = hello_and_goodbye.hello_output\n    }\n}",
    "options": "{\n\n}"
  },
  "calls": {
    "main_workflow.hello_and_goodbye": [{
      "executionStatus": "Done",
      "subWorkflowMetadata": {
        "workflowName": "hello_and_goodbye",
        "calls": {
          "sub.hello_and_goodbye.hello": [{
            "executionStatus": "Done",
            "stdout": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-hello/execution/stdout",
            "shardIndex": -1,
            "outputs": {
              "salutation": "Hello sub world!"
            },
            "runtimeAttributes": {
              "docker": "ubuntu:latest",
              "failOnStderr": false,
              "continueOnReturnCode": "0"
            },
            "cache": {
              "allowResultReuse": true
            },
            "Effective call caching mode": "CallCachingOff",
            "inputs": {
              "addressee": "sub world"
            },
            "returnCode": 0,
            "jobId": "49830",
            "backend": "Local",
            "end": "2016-11-17T14:13:40.712-05:00",
            "stderr": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-hello/execution/stderr",
            "callRoot": "cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-hello",
            "attempt": 1,
            "executionEvents": [{
              "startTime": "2016-11-17T14:13:39.240-05:00",
              "description": "Pending",
              "endTime": "2016-11-17T14:13:39.240-05:00"
            }, {
              "startTime": "2016-11-17T14:13:39.240-05:00",
              "description": "RequestingExecutionToken",
              "endTime": "2016-11-17T14:13:39.240-05:00"
            }, {
              "startTime": "2016-11-17T14:13:39.240-05:00",
              "description": "PreparingJob",
              "endTime": "2016-11-17T14:13:39.243-05:00"
            }, {
              "startTime": "2016-11-17T14:13:39.243-05:00",
              "description": "RunningJob",
              "endTime": "2016-11-17T14:13:40.704-05:00"
            }, {
              "startTime": "2016-11-17T14:13:40.704-05:00",
              "description": "UpdatingJobStore",
              "endTime": "2016-11-17T14:13:40.712-05:00"
            }],
            "start": "2016-11-17T14:13:39.239-05:00"
          }],
          "sub.hello_and_goodbye.goodbye": [{
            "executionStatus": "Done",
            "stdout": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-goodbye/execution/stdout",
            "shardIndex": -1,
            "outputs": {
              "salutation": "Goodbye sub world!"
            },
            "runtimeAttributes": {
              "docker": "ubuntu:latest",
              "failOnStderr": false,
              "continueOnReturnCode": "0"
            },
            "cache": {
              "allowResultReuse": true
            },
            "Effective call caching mode": "CallCachingOff",
            "inputs": {
              "addressee": "sub world"
            },
            "returnCode": 0,
            "jobId": "49831",
            "backend": "Local",
            "end": "2016-11-17T14:13:41.115-05:00",
            "stderr": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-goodbye/execution/stderr",
            "callRoot": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11/call-goodbye",
            "attempt": 1,
            "executionEvents": [{
              "startTime": "2016-11-17T14:13:39.240-05:00",
              "description": "Pending",
              "endTime": "2016-11-17T14:13:39.240-05:00"
            }, {
              "startTime": "2016-11-17T14:13:39.240-05:00",
              "description": "RequestingExecutionToken",
              "endTime": "2016-11-17T14:13:39.240-05:00"
            }, {
              "startTime": "2016-11-17T14:13:39.240-05:00",
              "description": "PreparingJob",
              "endTime": "2016-11-17T14:13:39.243-05:00"
            }, {
              "startTime": "2016-11-17T14:13:39.243-05:00",
              "description": "RunningJob",
              "endTime": "2016-11-17T14:13:41.112-05:00"
            }, {
              "startTime": "2016-11-17T14:13:41.112-05:00",
              "description": "UpdatingJobStore",
              "endTime": "2016-11-17T14:13:41.115-05:00"
            }],
            "start": "2016-11-17T14:13:39.239-05:00"
          }]
        },
        "outputs": {
          "goodbye_output": "Goodbye sub world!",
          "hello_output": "Hello sub world!"
        },
        "workflowRoot": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd/call-hello_and_goodbye/hello_and_goodbye/a6365f91-c807-465a-9186-a5d3da98fe11",
        "id": "a6365f91-c807-465a-9186-a5d3da98fe11",
        "inputs": {
          "hello_and_goodbye_input": "sub world"
        },
        "status": "Succeeded",
        "parentWorkflowId": "1d919bd4-d046-43b0-9918-9964509689dd",
        "end": "2016-11-17T14:13:41.116-05:00",
        "start": "2016-11-17T14:13:39.236-05:00"
      },
      "shardIndex": -1,
      "outputs": {
        "goodbye_output": "Goodbye sub world!",
        "hello_output": "Hello sub world!"
      },
      "inputs": {
        "hello_and_goodbye_input": "sub world"
      },
      "end": "2016-11-17T14:13:41.117-05:00",
      "attempt": 1,
      "start": "2016-11-17T14:13:39.236-05:00"
    }]
  },
  "outputs": {
    "main_output": "Hello sub world!"
  },
  "workflowRoot": "/cromwell-executions/main_workflow/1d919bd4-d046-43b0-9918-9964509689dd",
  "id": "1d919bd4-d046-43b0-9918-9964509689dd",
  "inputs": {

  },
  "submission": "2016-11-17T14:13:39.104-05:00",
  "status": "Succeeded",
  "end": "2016-11-17T14:13:41.120-05:00",
  "start": "2016-11-17T14:13:39.204-05:00"
}
```

# REST API

The `server` subcommand on the executable JAR will start an HTTP server which can accept WDL files to run as well as check status and output of existing workflows.

The following sub-sections define which HTTP Requests the web server can accept and what they will return.  Example HTTP requests are given in [HTTPie](https://github.com/jkbrzt/httpie) and [cURL](https://curl.haxx.se/)

## REST API Versions

All web server requests include an API version in the url. The current version is `v1`.

## POST /api/workflows/:version

This endpoint accepts a POST request with a `multipart/form-data` encoded body.  The form fields that may be included are:

* `wdlSource` - *Required* Contains the WDL file to submit for execution.
* `workflowInputs` - *Optional* JSON file containing the inputs.  A skeleton file can be generated from [wdltool](https://github.com/broadinstitute/wdltool) using the "inputs" subcommand.
* `workflowInputs_2` - *Optional* JSON file containing the inputs.
* `workflowInputs_3` - *Optional* JSON file containing the inputs.
* `workflowInputs_4` - *Optional* JSON file containing the inputs.
* `workflowInputs_5` - *Optional* JSON file containing the inputs.
* `workflowOptions` - *Optional* JSON file containing options for this workflow execution.  See the [run](#run) CLI sub-command for some more information about this.
* `wdlDependencies` - *Optional* ZIP file containing WDL files that are used to resolve import statements.

Regarding the workflowInputs parameter, in case of key conflicts between multiple input JSON files, higher values of x in workflowInputs_x override lower values. For example, an input specified in workflowInputs_3 will override an input with the same name in workflowInputs or workflowInputs_2.
Similarly, an input key specified in workflowInputs_5 will override an identical input key in any other input file.

Additionally, although Swagger has a limit of 5 JSON input files, the REST endpoint itself can accept an unlimited number of JSON input files.


cURL:

```
$ curl -v "localhost:8000/api/workflows/v1" -F wdlSource=@src/main/resources/3step.wdl -F workflowInputs=@test.json
```

HTTPie:

```
$ http --print=hbHB --form POST localhost:8000/api/workflows/v1 wdlSource=@src/main/resources/3step.wdl workflowInputs@inputs.json
```

Request:

```
POST /api/workflows/v1 HTTP/1.1
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

To specify workflow options as well:


cURL:

```
$ curl -v "localhost:8000/api/workflows/v1" -F wdlSource=@wdl/jes0.wdl -F workflowInputs=@wdl/jes0.json -F workflowOptions=@options.json
```

HTTPie:

```
http --print=HBhb --form POST http://localhost:8000/api/workflows/v1 wdlSource=@wdl/jes0.wdl workflowInputs@wdl/jes0.json workflowOptions@options.json
```

Request (some parts truncated for brevity):

```
POST /api/workflows/v1 HTTP/1.1
Accept: */*
Accept-Encoding: gzip, deflate
Connection: keep-alive
Content-Length: 1472
Content-Type: multipart/form-data; boundary=f3fd038395644de596c460257626edd7
Host: localhost:8000
User-Agent: HTTPie/0.9.2

--f3fd038395644de596c460257626edd7
Content-Disposition: form-data; name="wdlSource"

task x { ... }
task y { ... }
task z { ... }

workflow myworkflow {
  call x
  call y
  call z {
    input: example="gs://my-bucket/cromwell-executions/myworkflow/example.txt", int=3000
  }
}

--f3fd038395644de596c460257626edd7
Content-Disposition: form-data; name="workflowInputs"; filename="jes0.json"

{
  "myworkflow.x.x": "100"
}

--f3fd038395644de596c460257626edd7
Content-Disposition: form-data; name="workflowOptions"; filename="options.json"

{
  "jes_gcs_root": "gs://myworkflow-dev/workflows"
}

--f3fd038395644de596c460257626edd7--
```

## POST /api/workflows/:version/batch

This endpoint accepts a POST request with a `multipart/form-data`
encoded body.  The form fields that may be included are:

* `wdlSource` - *Required* Contains the WDL file to submit for
execution.
* `workflowInputs` - *Required* JSON file containing the inputs in a
JSON array. A skeleton file for a single inputs json element can be
generated from [wdltool](https://github.com/broadinstitute/wdltool)
using the "inputs" subcommand. The orderded endpoint responses will
contain one workflow submission response for each input, respectively.
* `workflowOptions` - *Optional* JSON file containing options for this
workflow execution.  See the [run](#run) CLI sub-command for some more
information about this.
* `wdlDependencies` - *Optional* ZIP file containing WDL files that are used to resolve import statements. Applied equally to all workflowInput sets.

cURL:

```
$ curl -v "localhost:8000/api/workflows/v1/batch" -F wdlSource=@src/main/resources/3step.wdl -F workflowInputs=@test_array.json
```

HTTPie:

```
$ http --print=hbHB --form POST localhost:8000/api/workflows/v1/batch wdlSource=@src/main/resources/3step.wdl workflowInputs@inputs_array.json
```

Request:

```
POST /api/workflows/v1/batch HTTP/1.1
Accept: */*
Accept-Encoding: gzip, deflate
Connection: keep-alive
Content-Length: 750
Content-Type: multipart/form-data; boundary=64128d499e9e4616adea7d281f695dcb
Host: localhost:8000
User-Agent: HTTPie/0.9.2

--64128d499e9e4616adea7d281f695dcb
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

--64128d499e9e4616adea7d281f695dcb
Content-Disposition: form-data; name="workflowInputs"; filename="inputs_array.json"

[
    {
        "three_step.cgrep.pattern": "..."
    },
    {
        "three_step.cgrep.pattern": "..."
    }
]

--64128d499e9e4616adea7d281f695dcb--
```

Response:

```
HTTP/1.1 201 Created
Content-Length: 96
Content-Type: application/json; charset=UTF-8
Date: Tue, 02 Jun 2015 18:06:28 GMT
Server: spray-can/1.3.3

[
    {
        "id": "69d1d92f-3895-4a7b-880a-82535e9a096e",
        "status": "Submitted"
    },
    {
        "id": "69d1d92f-3895-4a7b-880a-82535e9a096f",
        "status": "Submitted"
    }
]
```

To specify workflow options as well:


cURL:

```
$ curl -v "localhost:8000/api/workflows/v1/batch" -F wdlSource=@wdl/jes0.wdl -F workflowInputs=@wdl/jes0_array.json -F workflowOptions=@options.json
```

HTTPie:

```
http --print=HBhb --form POST http://localhost:8000/api/workflows/v1/batch wdlSource=@wdl/jes0.wdl workflowInputs@wdl/jes0_array.json workflowOptions@options.json
```

Request (some parts truncated for brevity):

```
POST /api/workflows/v1/batch HTTP/1.1
Accept: */*
Accept-Encoding: gzip, deflate
Connection: keep-alive
Content-Length: 1492
Content-Type: multipart/form-data; boundary=f3fd038395644de596c460257626edd8
Host: localhost:8000
User-Agent: HTTPie/0.9.2

--f3fd038395644de596c460257626edd8
Content-Disposition: form-data; name="wdlSource"

task x { ... }
task y { ... }
task z { ... }

workflow myworkflow {
  call x
  call y
  call z {
    input: example="gs://my-bucket/cromwell-executions/myworkflow/example.txt", int=3000
  }
}

--f3fd038395644de596c460257626edd8
Content-Disposition: form-data; name="workflowInputs"; filename="jes0_array.json"

[
  {
    "myworkflow.x.x": "100"
  }, {
    "myworkflow.x.x": "101"
  }
]

--f3fd038395644de596c460257626edd8
Content-Disposition: form-data; name="workflowOptions"; filename="options.json"

{
  "jes_gcs_root": "gs://myworkflow-dev/workflows"
}

--f3fd038395644de596c460257626edd8--
```


## GET /api/workflows/:version/query

This endpoint allows for querying workflows based on the following criteria:

* `name`
* `id`
* `status`
* `start` (start datetime with mandatory offset)
* `end` (end datetime with mandatory offset)
* `page` (page of results)
* `pagesize` (# of results per page)

Names, ids, and statuses can be given multiple times to include
workflows with any of the specified names, ids, or statuses. When
multiple names are specified, any workflow matching one of the names
will be returned. The same is true for multiple ids or statuses. When
different types of criteria are specified, for example names and
statuses, the results must match both the one of the specified names and
one of the statuses. Using page and pagesize will enable server side pagination.

Valid statuses are `Submitted`, `Running`, `Aborting`, `Aborted`, `Failed`, and `Succeeded`.  `start` and `end` should
be in [ISO8601 datetime](https://en.wikipedia.org/wiki/ISO_8601) format with *mandatory offset* and `start` cannot be after `end`.

cURL:

```
$ curl "http://localhost:8000/api/workflows/v1/query?start=2015-11-01T00%3A00%3A00-04%3A00&end=2015-11-04T00%3A00%3A00-04%3A00&status=Failed&status=Succeeded&page=1&pagesize=10"
```

HTTPie:

```
$ http "http://localhost:8000/api/workflows/v1/query?start=2015-11-01T00%3A00%3A00-04%3A00&end=2015-11-04T00%3A00%3A00-04%3A00&status=Failed&status=Succeeded&page=1&pagesize=10"
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 133
Content-Type: application/json; charset=UTF-8
Date: Tue, 02 Jun 2015 18:06:56 GMT
Server: spray-can/1.3.3

{
  "results": [
    {
      "name": "w",
      "id": "fdfa8482-e870-4528-b639-73514b0469b2",
      "status": "Succeeded",
      "end": "2015-11-01T07:45:52.000-05:00",
      "start": "2015-11-01T07:38:57.000-05:00"
    },
    {
      "name": "hello",
      "id": "e69895b1-42ed-40e1-b42d-888532c49a0f",
      "status": "Succeeded",
      "end": "2015-11-01T07:45:30.000-05:00",
      "start": "2015-11-01T07:38:58.000-05:00"
    },
    {
      "name": "crasher",
      "id": "ed44cce4-d21b-4c42-b76d-9d145e4d3607",
      "status": "Failed",
      "end": "2015-11-01T07:45:44.000-05:00",
      "start": "2015-11-01T07:38:59.000-05:00"
    }
  ]
}
```

Query data is refreshed from raw data periodically according to the configuration value `services.MetadataService.metadata-summary-refresh-interval`.
This interval represents the duration between the end of one summary refresh sweep and the beginning of the next sweep.  If not specified the
refresh interval will default to 2 seconds.  To turn off metadata summary refresh, specify an infinite refresh interval value with "Inf".

```
services {
  MetadataService {
    metadata-summary-refresh-interval = "10 seconds"
  }
}
```

## POST /api/workflows/:version/query

This endpoint allows for querying workflows based on the same criteria
as [GET /api/workflows/:version/query](#get-apiworkflowsversionquery).

Instead of specifying query parameters in the URL, the parameters
must be sent via the POST body. The request content type must be
`application/json`. The json should be a list of objects. Each json
object should contain a different criterion.

cURL:

```
$ curl -X POST --header "Content-Type: application/json" -d "[{\"start\": \"2015-11-01T00:00:00-04:00\"}, {\"end\": \"2015-11-04T00:00:00-04:00\"}, {\"status\": \"Failed\"}, {\"status\": \"Succeeded\"}, {\"page\": \"1\"}, {\"pagesize\": \"10\"}]" "http://localhost:8000/api/workflows/v1/query"
```

HTTPie:

```
$ echo "[{\"start\": \"2015-11-01T00:00:00-04:00\"}, {\"end\": \"2015-11-04T00:00:00-04:00\"}, {\"status\": \"Failed\"}, {\"status\": \"Succeeded\"}, {\"page\": \"1\"}, {\"pagesize\": \"10\"}]" | http "http://localhost:8000/api/workflows/v1/query"
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 133
Content-Type: application/json; charset=UTF-8
Date: Tue, 02 Jun 2015 18:06:56 GMT
Server: spray-can/1.3.3

{
  "results": [
    {
      "name": "w",
      "id": "fdfa8482-e870-4528-b639-73514b0469b2",
      "status": "Succeeded",
      "end": "2015-11-01T07:45:52.000-05:00",
      "start": "2015-11-01T07:38:57.000-05:00"
    },
    {
      "name": "hello",
      "id": "e69895b1-42ed-40e1-b42d-888532c49a0f",
      "status": "Succeeded",
      "end": "2015-11-01T07:45:30.000-05:00",
      "start": "2015-11-01T07:38:58.000-05:00"
    },
    {
      "name": "crasher",
      "id": "ed44cce4-d21b-4c42-b76d-9d145e4d3607",
      "status": "Failed",
      "end": "2015-11-01T07:45:44.000-05:00",
      "start": "2015-11-01T07:38:59.000-05:00"
    }
  ]
}
```

## GET /api/workflows/:version/:id/status

cURL:

```
$ curl http://localhost:8000/api/workflows/v1/69d1d92f-3895-4a7b-880a-82535e9a096e/status
```

HTTPie:

```
$ http http://localhost:8000/api/workflows/v1/69d1d92f-3895-4a7b-880a-82535e9a096e/status
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

## GET /api/workflows/:version/:id/outputs

cURL:

```
$ curl http://localhost:8000/api/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/outputs
```

HTTPie:

```
$ http http://localhost:8000/api/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/outputs
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

## GET /api/workflows/:version/:id/timing

This endpoint is meant to be used in a web browser.  It will show a Gantt Chart of a particular workflow.  The bars in the chart represent start and end times for individual task invocations.

![Timing diagram](http://i.imgur.com/EOE2HoL.png)

## GET /api/workflows/:version/:id/logs

This will return paths to the standard out and standard error files that were generated during the execution of all calls in a workflow.

A call has one or more standard out and standard error logs, depending on if the call was scattered or not. In the latter case, one log is provided for each instance of the call that has been run.

cURL:

```
$ curl http://localhost:8000/api/workflows/v1/b3e45584-9450-4e73-9523-fc3ccf749848/logs
```

HTTPie:

```
$ http http://localhost:8000/api/workflows/v1/b3e45584-9450-4e73-9523-fc3ccf749848/logs
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
        "call.ps": [
            {
                "stderr": "/home/user/test/b3e45584-9450-4e73-9523-fc3ccf749848/call-ps/stderr6126967977036995110.tmp",
                "stdout": "/home/user/test/b3e45584-9450-4e73-9523-fc3ccf749848/call-ps/stdout6128485235785447571.tmp"
            }
        ],
        "call.cgrep": [
            {
                "stderr": "/home/user/test/b3e45584-9450-4e73-9523-fc3ccf749848/call-cgrep/stderr6126967977036995110.tmp",
                "stdout": "/home/user/test/b3e45584-9450-4e73-9523-fc3ccf749848/call-cgrep/stdout6128485235785447571.tmp"
            }
        ],
        "call.wc": [
            {
                "stderr": "/home/user/test/b3e45584-9450-4e73-9523-fc3ccf749848/call-wc/stderr6126967977036995110.tmp",
                "stdout": "/home/user/test/b3e45584-9450-4e73-9523-fc3ccf749848/call-wc/stdout6128485235785447571.tmp"
            }
        ]
    }
}
```

## GET /api/workflows/:version/:id/metadata

This endpoint returns a superset of the data from #get-workflowsversionidlogs in essentially the same format
(i.e. shards are accounted for by an array of maps, in the same order as the shards). 
In addition to shards, every attempt that was made for this call will have its own object as well, in the same order as the attempts.
Workflow metadata includes submission, start, and end datetimes, as well as status, inputs and outputs.
Call-level metadata includes inputs, outputs, start and end datetime, backend-specific job id,
return code, stdout and stderr.  Date formats are ISO with milliseconds.

Accepted parameters are:

* `includeKey` Optional repeated string value, specifies what metadata
  keys to include in the output, matched as a prefix string. Keys that
  are not specified are filtered out. The call keys `attempt` and
  `shardIndex` will always be included. May not be used with
  `excludeKey`.

* `excludeKey` Optional repeated string value, specifies what metadata
  keys to exclude from the output, matched as a prefix string. Keys that
  are specified are filtered out. The call keys `attempt` and
  `shardIndex` will always be included. May not be used with
  `includeKey`.

cURL:

```
$ curl http://localhost:8000/api/workflows/v1/b3e45584-9450-4e73-9523-fc3ccf749848/metadata
```

HTTPie:

```
$ http http://localhost:8000/api/workflows/v1/b3e45584-9450-4e73-9523-fc3ccf749848/metadata
```

Response:
```
HTTP/1.1 200 OK
Server spray-can/1.3.3 is not blacklisted
Server: spray-can/1.3.3
Date: Thu, 01 Oct 2015 22:18:07 GMT
Content-Type: application/json; charset=UTF-8
Content-Length: 7286
{
  "workflowName": "sc_test",
  "calls": {
    "sc_test.do_prepare": [
      {
        "executionStatus": "Done",
        "stdout": "/home/jdoe/cromwell/cromwell-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_prepare/stdout",
        "shardIndex": -1,
        "outputs": {
          "split_files": [
            "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_prepare/temp_aa",
            "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_prepare/temp_ad"
          ]
        },
        "inputs": {
          "input_file": "/home/jdoe/cromwell/11.txt"
        },
        "runtimeAttributes": {
            "failOnStderr": "true",
            "continueOnReturnCode": "0"
        },
        "returnCode": 0,
        "backend": "Local",
        "end": "2016-02-04T13:47:56.000-05:00",
        "stderr": "/home/jdoe/cromwell/cromwell-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_prepare/stderr",
        "attempt": 1,
        "executionEvents": [],
        "start": "2016-02-04T13:47:55.000-05:00"
      }
    ],
    "sc_test.do_scatter": [
      {
        "executionStatus": "Preempted",
        "stdout": "/home/jdoe/cromwell/cromwell-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-0/stdout",
        "shardIndex": 0,
        "outputs": {},
        "runtimeAttributes": {
           "failOnStderr": "true",
           "continueOnReturnCode": "0"
        },
        "inputs": {
          "input_file": "f"
        },
        "backend": "Local",
        "end": "2016-02-04T13:47:56.000-05:00",
        "stderr": "/home/jdoe/cromwell/cromwell-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-0/stderr",
        "attempt": 1,
        "executionEvents": [],
        "start": "2016-02-04T13:47:56.000-05:00"
      },
      {
        "executionStatus": "Done",
        "stdout": "/home/jdoe/cromwell/cromwell-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-0/attempt-2/stdout",
        "shardIndex": 0,
        "outputs": {
          "count_file": "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-0/attempt-2/output.txt"
        },
        "runtimeAttributes": {
           "failOnStderr": "true",
           "continueOnReturnCode": "0"
        },
        "inputs": {
          "input_file": "f"
        },
        "returnCode": 0,
        "end": "2016-02-04T13:47:56.000-05:00",
        "stderr": "/home/jdoe/cromwell/cromwell-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-0/attempt-2/stderr",
        "attempt": 2,
        "executionEvents": [],
        "start": "2016-02-04T13:47:56.000-05:00"
      },
      {
        "executionStatus": "Done",
        "stdout": "/home/jdoe/cromwell/cromwell-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-1/stdout",
        "shardIndex": 1,
        "outputs": {
          "count_file": "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-1/output.txt"
        },
        "runtimeAttributes": {
           "failOnStderr": "true",
           "continueOnReturnCode": "0"
        },
        "inputs": {
          "input_file": "f"
        },
        "returnCode": 0,
        "backend": "Local",
        "end": "2016-02-04T13:47:56.000-05:00",
        "stderr": "/home/jdoe/cromwell/cromwell-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-1/stderr",
        "attempt": 1,
        "executionEvents": [],
        "start": "2016-02-04T13:47:56.000-05:00"
      }
    ],
    "sc_test.do_gather": [
      {
        "executionStatus": "Done",
        "stdout": "/home/jdoe/cromwell/cromwell-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_gather/stdout",
        "shardIndex": -1,
        "outputs": {
          "sum": 12
        },
        "runtimeAttributes": {
           "failOnStderr": "true",
           "continueOnReturnCode": "0"
        },
        "inputs": {
          "input_files": [
            "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-0/attempt-2/output.txt",
            "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-1/output.txt"
          ]
        },
        "returnCode": 0,
        "backend": "Local",
        "end": "2016-02-04T13:47:57.000-05:00",
        "stderr": "/home/jdoe/cromwell/cromwell-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_gather/stderr",
        "attempt": 1,
        "executionEvents": [],
        "start": "2016-02-04T13:47:56.000-05:00"
      }
    ]
  },
  "outputs": {
    "sc_test.do_gather.sum": 12,
    "sc_test.do_prepare.split_files": [
      "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_prepare/temp_aa",
      "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_prepare/temp_ad"
    ],
    "sc_test.do_scatter.count_file": [
      "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-0/attempt-2/output.txt",
      "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-1/output.txt"
    ]
  },
  "id": "8e592ed8-ebe5-4be0-8dcb-4073a41fe180",
  "inputs": {
    "sc_test.do_prepare.input_file": "/home/jdoe/cromwell/11.txt"
  },
  "submission": "2016-02-04T13:47:55.000-05:00",
  "status": "Succeeded",
  "end": "2016-02-04T13:47:57.000-05:00",
  "start": "2016-02-04T13:47:55.000-05:00"
}
```

cURL:

```
$ curl "http://localhost:8000/api/workflows/v1/b3e45584-9450-4e73-9523-fc3ccf749848/metadata?includeKey=inputs&includeKey=outputs"
```

HTTPie:

```
$ http "http://localhost:8000/api/workflows/v1/b3e45584-9450-4e73-9523-fc3ccf749848/metadata?includeKey=inputs&includeKey=outputs"
```

Response:
```
HTTP/1.1 200 OK
Server spray-can/1.3.3 is not blacklisted
Server: spray-can/1.3.3
Date: Thu, 01 Oct 2015 22:19:07 GMT
Content-Type: application/json; charset=UTF-8
Content-Length: 4286
{
  "calls": {
    "sc_test.do_prepare": [
      {
        "shardIndex": -1,
        "outputs": {
          "split_files": [
            "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_prepare/temp_aa",
            "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_prepare/temp_ad"
          ]
        },
        "inputs": {
          "input_file": "/home/jdoe/cromwell/11.txt"
        },
        "attempt": 1
      }
    ],
    "sc_test.do_scatter": [
      {
        "shardIndex": 0,
        "outputs": {},
        "inputs": {
          "input_file": "f"
        },
        "attempt": 1
      },
      {
        "shardIndex": 0,
        "outputs": {
          "count_file": "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-0/attempt-2/output.txt"
        },
        "inputs": {
          "input_file": "f"
        },
        "attempt": 2
      },
      {
        "shardIndex": 1,
        "outputs": {
          "count_file": "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-1/output.txt"
        },
        "inputs": {
          "input_file": "f"
        },
        "attempt": 1
      }
    ],
    "sc_test.do_gather": [
      {
        "shardIndex": -1,
        "outputs": {
          "sum": 12
        }
        "inputs": {
          "input_files": [
            "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-0/attempt-2/output.txt",
            "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-1/output.txt"
          ]
        },
        "attempt": 1
      }
    ]
  },
  "outputs": {
    "sc_test.do_gather.sum": 12,
    "sc_test.do_prepare.split_files": [
      "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_prepare/temp_aa",
      "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_prepare/temp_ad"
    ],
    "sc_test.do_scatter.count_file": [
      "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-0/attempt-2/output.txt",
      "/home/jdoe/cromwell/cromwell-test-executions/sc_test/8e592ed8-ebe5-4be0-8dcb-4073a41fe180/call-do_scatter/shard-1/output.txt"
    ]
  },
  "id": "8e592ed8-ebe5-4be0-8dcb-4073a41fe180",
  "inputs": {
    "sc_test.do_prepare.input_file": "/home/jdoe/cromwell/11.txt"
  }
}
```

The `call` and `workflow` may optionally contain failures shaped like this:
```
"failures": [
  {
    "failure": "The failure message",
    "timestamp": "2016-02-25T10:49:02.066-05:00"
  }
]
```

## POST /api/workflows/:version/:id/abort

cURL:

```
$ curl -X POST http://localhost:8000/api/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/abort
```

HTTPie:

```
$ http POST http://localhost:8000/api/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/abort
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

## GET /api/workflows/:version/backends

This endpoint returns a list of the backends supported by the server as well as the default backend.

cURL:
```
$ curl http://localhost:8000/api/workflows/v1/backends
```

HTTPie:
```
$ http http://localhost:8000/api/workflows/v1/backends
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 379
Content-Type: application/json; charset=UTF-8
Date: Mon, 03 Aug 2015 17:11:28 GMT
Server: spray-can/1.3.3

{
  "supportedBackends": ["JES", "LSF", "Local", "SGE"],
  "defaultBackend": "Local"
}
```

## GET /api/engine/:version/stats

This endpoint returns some basic statistics on the current state of the engine. At the moment that includes the number of running workflows and the number of active jobs. 

cURL:
```
$ curl http://localhost:8000/api/engine/v1/stats
```

HTTPie:
```
$ http http://localhost:8000/api/engine/v1/stats
```

Response:
```
"date": "Sun, 18 Sep 2016 14:38:11 GMT",
"server": "spray-can/1.3.3",
"content-length": "33",
"content-type": "application/json; charset=UTF-8"

{
  "workflows": 3,
  "jobs": 10
}
```

## GET /api/engine/:version/version

This endpoint returns the version of the Cromwell engine.

cURL:
```
$ curl http://localhost:8000/api/engine/v1/version
```

HTTPie:
```
$ http http://localhost:8000/api/engine/v1/version
```

Response:
```
"date": "Sun, 18 Sep 2016 14:38:11 GMT",
"server": "spray-can/1.3.3",
"content-length": "33",
"content-type": "application/json; charset=UTF-8"

{
  "cromwell": 23-8be799a-SNAP
}
```



## Error handling
Requests that Cromwell can't process return a failure in the form of a JSON response respecting the following JSON schema:

```
{
  "$schema": "http://json-schema.org/draft-04/schema#",
  "description": "Error response schema",
  "type": "object",
  "properties": {
    "status": {
      "enum": [ "fail", "error"]
    },
    "message": {
      "type": "string"
    },
    "errors": {
      "type": "array",
      "minItems": 1,
      "items": { "type": "string" },
      "uniqueItems": true
    }
  },
  "required": ["status", "message"]
}
```

The `status` field can take two values:
> "fail" means that the request was invalid and/or data validation failed. "fail" status is most likely returned with a 4xx HTTP Status code.
e.g.

```
{
  "status": "fail",
  "message": "Workflow input processing failed.",
  "errors": [
    "Required workflow input 'helloworld.input' not specified."
  ]
}
```

> "error" means that an error occurred while processing the request. "error" status is most likely returned with a 5xx HTTP Status code.
e.g.

```
{
  "status": "error",
  "message": "Connection to the database failed."
}
```

The `message` field contains a short description of the error.

The `errors` field is optional and may contain additional information about why the request failed.

# Developer

## Generating table of contents on Markdown files

```
$ pip install mdtoc
$ mdtoc --check-links README.md
```

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
