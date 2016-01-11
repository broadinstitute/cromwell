[![Build Status](https://travis-ci.org/broadinstitute/cromwell.svg?branch=develop)](https://travis-ci.org/broadinstitute/cromwell?branch=develop)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/cromwell/badge.svg?branch=develop)](https://coveralls.io/r/broadinstitute/cromwell?branch=develop)

Cromwell
========

A [Workflow Management System](https://en.wikipedia.org/wiki/Workflow_management_system) geared towards [scientific workflows](https://en.wikipedia.org/wiki/Scientific_workflow_system). Cromwell is open sourced under the BSD 3-Clause license.

<!---toc start-->

* [Mailing List](#mailing-list)
* [Requirements](#requirements)
* [Building](#building)
* [Installing](#installing)
* [Command Line Usage](#command-line-usage)
  * [validate](#validate)
  * [inputs](#inputs)
  * [run](#run)
  * [parse](#parse)
  * [highlight](#highlight)
  * [server](#server)
* [Getting Started with WDL](#getting-started-with-wdl)
* [Configuring Cromwell](#configuring-cromwell)
  * [Database](#database)
* [Backends](#backends)
  * [Backend Filesystems](#backend-filesystems)
    * [Shared Local Filesystem](#shared-local-filesystem)
    * [Google Cloud Storage Filesystem](#google-cloud-storage-filesystem)
  * [Local Backend](#local-backend)
  * [Sun GridEngine Backend](#sun-gridengine-backend)
  * [Google JES Backend](#google-jes-backend)
    * [Data Localization](#data-localization)
    * [Docker](#docker)
    * [Monitoring](#monitoring)
* [Runtime Attributes](#runtime-attributes)
  * [continueOnReturnCode](#continueonreturncode)
  * [cpu](#cpu)
  * [defaultDisks](#defaultdisks)
  * [defaultZones](#defaultzones)
  * [docker](#docker)
  * [failOnStderr](#failonstderr)
  * [memory](#memory)
  * [preemptible](#preemptible)
* [Logging](#logging)
* [Workflow Options](#workflow-options)
* [Call Caching](#call-caching)
* [REST API](#rest-api)
  * [REST API Versions](#rest-api-versions)
  * [POST /api/workflows/:version](#post-apiworkflowsversion)
  * [GET /api/workflows/:version/query](#get-apiworkflowsversionquery)
  * [GET /api/workflows/:version/:id/status](#get-apiworkflowsversionidstatus)
  * [GET /api/workflows/:version/:id/outputs](#get-apiworkflowsversionidoutputs)
  * [GET /api/workflows/:version/:id/timing](#get-apiworkflowsversionidtiming)
  * [GET /api/workflows/:version/:id/outputs/:call](#get-apiworkflowsversionidoutputscall)
  * [GET /api/workflows/:version/:id/logs/:call](#get-apiworkflowsversionidlogscall)
  * [GET /api/workflows/:version/:id/logs](#get-apiworkflowsversionidlogs)
  * [GET /api/workflows/:version/:id/metadata](#get-apiworkflowsversionidmetadata)
  * [POST /api/workflows/:version/:id/abort](#post-apiworkflowsversionidabort)
  * [POST /api/workflows/:version/:id/call-caching](#post-apiworkflowsversionidcall-caching)
  * [POST /api/workflows/:version/:id/call-caching/:call](#post-apiworkflowsversionidcall-cachingcall)
* [Developer](#developer)
  * [Generating table of contents on Markdown files](#generating-table-of-contents-on-markdown-files)
  * [Generating and Hosting ScalaDoc](#generating-and-hosting-scaladoc)

<!---toc end-->

# Mailing List

The [Cromwell Mailing List](https://groups.google.com/a/broadinstitute.org/forum/?hl=en#!forum/cromwell) is cromwell@broadinstitute.org.

If you have any questions, suggestions or support issues please send them to this list. To subscribe you can either join via the link above or send an email to cromwell+subscribe@broadinstitute.org.

# Requirements

The following is the toolchain used for development of Cromwell.  Other versions may work, but these are recommended.

* [Scala 2.11.7](http://www.scala-lang.org/news/2.11.7)
* [SBT 0.13.8](https://github.com/sbt/sbt/releases/tag/v0.13.8)
* [Java 8](http://www.oracle.com/technetwork/java/javase/overview/java8-2100321.html)

# Building

`sbt assembly` will build a runnable JAR in `target/scala-2.11/`

Tests are run via `sbt test`.  Note that the tests do require Docker to be running.  To test this out while downloading the Ubuntu image that is required for tests, run `docker pull ubuntu:latest` prior to running `sbt test`

# Installing

OS X users can install Cromwell with Homebrew: `brew install cromwell`.

# Command Line Usage

Run the JAR file with no arguments to get the usage message:

```
$ java -jar cromwell.jar
java -jar cromwell.jar <action> <parameters>

Actions:
run <WDL file> [<JSON inputs file> [<JSON workflow options>
  [<OUTPUT workflow metadata>]]]

  Given a WDL file and JSON file containing the value of the
  workflow inputs, this will run the workflow locally and
  print out the outputs in JSON format.  The workflow
  options file specifies some runtime configuration for the
  workflow (see README for details).  The workflow metadata
  output is an optional file path to output the metadata.
  Use a single dash ("-") to skip optional files. Ex:
    run noinputs.wdl - - metadata.json

server

  Starts a web server on port 8000.  See the web server
  documentation for more details about the API endpoints.

parse <WDL file>

  Compares a WDL file against the grammar and prints out an
  abstract syntax tree if it is valid, and a syntax error
  otherwise.  Note that higher-level AST checks are not done
  via this sub-command and the 'validate' subcommand should
  be used for full validation

validate <WDL file>

  Performs full validation of the WDL file including syntax
  and semantic checking

inputs <WDL file>

  Print a JSON skeleton file of the inputs needed for this
  workflow.  Fill in the values in this JSON document and
  pass it in to the 'run' subcommand.

highlight <WDL file> <html|console>

  Reformats and colorizes/tags a WDL file. The second
  parameter is the output type.  "html" will output the WDL
  file with <span> tags around elements.  "console" mode
  will output colorized text to the terminal
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

Only a few workflow options are available currently and are all to be used with the JES backend. See the section on the [JES backend](#google-jes-backend) for more details.

```
$ java -jar cromwell.jar run my_jes_wf.wdl my_jes_wf.json wf_options.json
```

The fourth, optional parameter to the 'run' subcommand is a path where the workflow metadata will be written.  By default, no workflow metadata will be written.

```
$ java -jar cromwell.jar run my_wf.wdl - - my_wf.metadata.json
... play-by-play output ...
$ cat my_wf.metadata.json
{
  "calls": {
    "example.my_task": [{
      "executionStatus": "Done",
      "stdout": "/Users/cromwell/cromwell-executions/example/22b6f829-e2f9-4813-9d20-3328669c786b/call-my_task/stdout",
      "outputs": {
        "result": "my example output"
      },
      "inputs": {
      },
      "returnCode": 0,
      "backend": "Local",
      "end": "2015-10-29T03:16:51.732-03:00",
      "stderr": "/Users/cromwell/cromwell-executions/example/22b6f829-e2f9-4813-9d20-3328669c786b/call-my_task/stderr",
      "start": "2015-10-29T03:16:51.213-03:00"
      "executionEvents": [{
        "description": "running docker",
        "startTime": "2015-10-29T03:16:51.213-03:00",
        "endTime": "2015-10-29T03:16:51.732-03:00"
      }]
    }]
  },
  "outputs": {
    "example.my_task.result": "my = /root/22b6f829-e2f9-4813-9d20-3328669c786b/call-my_task"
  },
  "id": "22b6f829-e2f9-4813-9d20-3328669c786b",
  "inputs": {
  },
  "submission": "2015-10-29T03:16:51.125-03:00",
  "status": "Succeeded",
  "end": "2015-10-29T03:16:51.740-03:00",
  "start": "2015-10-29T03:16:51.125-03:00"
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

For many examples on how to use WDL see [the WDL site](https://github.com/broadinstitute/wdl/tree/develop#getting-started-with-wdl)

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

## Database

Cromwell uses either an in-memory or MySQL database to track the execution of workflows and store outputs of task invocations.

By default, Cromwell uses an in-memory database which will only live for the duration of the JVM.  This provides a quick way to run workflows locally without having to set up MySQL, though it also makes workflow executions somewhat transient.

To configure Cromwell to instead point to a MySQL database, first create the empty database.  In the example below, the database name is `cromwell`.

Then, edit the configuration file `database` stanza, as follows:

```
database {
  config = main.mysql

  main {
    hsqldb {
      db.url = "jdbc:hsqldb:mem:${slick.uniqueSchema};shutdown=false;hsqldb.tx=mvcc"
      db.driver = "org.hsqldb.jdbcDriver"
      driver = "slick.driver.HsqldbDriver$"
      slick.createSchema = true
    }
    mysql {
      db.url = "jdbc:mysql://localhost:3306/cromwell"
      db.user = "root"
      db.password = ""
      db.driver = "com.mysql.jdbc.Driver"
      driver = "slick.driver.MySQLDriver$"
    }
  }

  test {
    ...
  }
}
```

To initially populate the tables, use the [Java MySQL Connector](https://dev.mysql.com/downloads/connector/j/) JAR file with [Liquibase](http://www.liquibase.org/) (installable via `brew install liquibase`):

```
liquibase --driver=com.mysql.jdbc.Driver \
          --classpath=${HOME}/.ivy2/cache/mysql/mysql-connector-java/jars/mysql-connector-java-5.1.35.jar \
          --changeLogFile=src/main/migrations/changelog.xml \
          --url="jdbc:mysql://localhost/cromwell" \
          --username="root" \
          --password="" \
          migrate
```

# Backends

A backend represents a way to run the user's command specified in the `task` section.  Currently three backends are supported:

* Local - Run jobs as subprocesses.  Supports launching in Docker containers.
* Sun GridEngine - Use `qsub` and job monitoring to run scripts.
* Google JES - Launch jobs on Google Compute Engine through the Job Execution Service (JES).

Backends are specified via the configuration option `backend.backend` which can accept the values: `sge`, `local`, and `jes` (e.g. `java -Dbackend.backend=sge`).

## Backend Filesystems

Each backend will utilize one filesystem to store the directory structure of an executed workflow.  Currently, the backend and the type of filesystem that the backend uses are tightly coupled.  In future versions of Cromwell, they may be more loosely coupled.

The backend/filesystem pairings are as follows:

* [Local Backend](#local-backend) uses the [Shared Local Filesystem](#shared-local-filesystem)
* [SGE Backend](#sun-gridengine-backend) uses the [Shared Local Filesystem](#shared-local-filesystem)
* [JES Backend](#google-jes-backend) uses the [Google Cloud Storage Filesystem](#google-cloud-storage-filesystem)

### Shared Local Filesystem

For the [local](#local-backend) and [Sun GridEngine](#sun-gridengine-backend) backends, the following is required of the underlying filesystem:

* (`local` backend) Subprocesses that Cromwell launches can use child directories that Cromwell creates as their CWD.  The subprocess must have write access to the directory that Cromwell assigns as its current working directory.
* (`sge` backend) Jobs launched with `qsub` can use directories that Cromwell creates as the working directory of the job, and write files to those directories.

The root directory that Cromwell uses for all workflows (`cromwell-root`) defaults to `./cromwell-executions`.  However, this is can be overwritten with the `-Dbackend.shared-filesystem.root=/your/path` option on the command line, or via [Cromwell's configuration file](#configuring-cromwell)

When cromwell runs a workflow, it first creates a directory `<cromwell-root>/<workflow_uuid>`.  This is called the `workflow_root` and it is the root directory for all activity in this workflow.

Each `call` has its own subdirectory located at `<workflow_root>/call-<call_name>`.  This is the `<call_dir>`.  Within this directory are special files written by the backend and they're supposed to be backend specific things though there are commonalities.  For example, having a `stdout` and `stderr` file is common among both backends and they both write a shell script file to the `<call_dir>` as well.  See the descriptions below for details about backend-specific files that are written to these directories.

An example of a workflow output directory would look like this:

```
cromwell-executions/
└── three_step
    └── df6363df-d812-4088-acfd-5b00ef3f5dcc
        ├── call-cgrep
        │   ├── cromwell_root
        │   │   └── cromwell
        │   │       └── cromwell-executions
        │   │           └── three_step
        │   │               └── df6363df-d812-4088-acfd-5b00ef3f5dcc
        │   │                   └── call-ps
        │   │                       └── stdout
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
            ├── cromwell_root
            │   └── cromwell
            │       └── cromwell-executions
            │           └── three_step
            │               └── df6363df-d812-4088-acfd-5b00ef3f5dcc
            │                   └── call-ps
            │                       └── stdout
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

This workflow output directory would be the result of running the above WDL file with Cromwell from the directory `/cromwell_root`.

In the above directory structure, you'll notice that the `call-cgrep` and `call-wc` sub-directories both contain a directory structure to point to the `stdout` file from the invocation of `ps`.  In these cases, that `stdout` file is a localized version of the one within `call-ps/stdout`.  By default both of those `stdout` files would be hard-links but they could also be symbolic links or copies of the file, depending on how Cromwell is configured (see below).  The directory structure is nested so deeply to avoid collisions.  For example, if either of these call invocations referenced two files called `stdout`, they'd collide if they were put into the same directory so the full directory structure is maintained.

Any input files to a call need to be localized into the `<call_dir>`.  There are a few localization strategies that Cromwell will try until one works.  Below is the default order specified in `application.conf` but this order can be overridden:

* `hard-link` - This will create a hard link (not symbolic) link to the file
* `soft-link` - Create a symbolic link to the file.  This strategy is not applicable for tasks which specify a Docker image and will be ignored.
* `copy` - Make a copy the file

These options can be overridden with command line options to Java.  For instance, to use the strategies `copy` and `hard-link`, in that order:

```
java -Dbackend.shared-filesystem.localization.0=copy -Dbackend.shared-filesystem.localization.1=hard-link cromwell.jar ...
```

Backends that use the shared filesystem can accept the following values for `File` variables:

* Local file system paths, either relative or absolute.  Relative paths are interpreted as relative to the current working directory of the Cromwell process.
* [Google Cloud Storage](https://cloud.google.com/storage/) URIs (e.g. `gs://my-bucket/x/y/z.txt`).  Any GCS URI will be downloaded locally

### Google Cloud Storage Filesystem

The Google Cloud Storage (GCS) Filesystem is only used for when a workflow is run on Google JES.  It uses the same directory structure as the [Shared Local Filesystem](#shared-local-filesystem), however it is rooted at one of the following GCS locations:

* If the `jes_gcs_root` [workflow option](#workflow-options) is set, this is used first.
* Otherwise, `backend.jes.baseExecutionBucket` in the [configuration file](#configuring-cromwell), which can also be set via `java -Dbackend.jes.baseExecutionBucket="gs://my-bucket/"`, will be used instead.

Google Cloud Storage URIs are the only acceptable values for `File` inputs for workflows using the JES backend.

## Local Backend

The local backend will simply launch a subprocess for each task invocation and wait for it to exit.

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

> **NOTE**: If you are using the local backend with Docker and Docker Machine on Mac OS X, by default Cromwell can only
> run from in any path under your home directory.
>
> The `-v` flag will only work if `<local_workflow_dir>` is within your home directory because VirtualBox with
> Docker Machine only exposes the home directory by default.  Any local path used in `-v` that is not within the user's
> home directory will silently be interpreted as references to paths on the VirtualBox VM.  This can manifest in
> Cromwell as tasks failing for odd reasons (like missing RC file)
>
> See https://docs.docker.com/engine/userguide/dockervolumes/ for more information on volume mounting in Docker.

## Sun GridEngine Backend

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

## Google JES Backend

Google JES (Job Execution Service) is a Docker-as-a-service from Google. JES has some [configuration](#configuring-cromwell) that needs to be set before it can be run.  Edit `src/main/resources/application.conf` and fill out the 'jes' stanza, e.g.

```hocon
backend {
  backend = "jes"

  jes {
    // Google project
    project = "broad-dsde-dev"

    // Location to store workflow results, must be a gs:// URL
    baseExecutionBucket = "gs://your-bucket/cromwell-executions"

    // Root URL for the API
    endpointUrl = "https://genomics.googleapis.com/"

    // Polling for completion backs-off gradually for slower-running jobs.
    // This is the maximum polling interval (in seconds):
    maximumPollingInterval = 600
  }

  ...
}
```

It is also necessary to fill out the `google` stanza in the configuration file. This stanza will set up the service / user account that Cromwell uses to write certain files to GCS as well as run jobs.

```hocon
google {
  applicationName = "cromwell"
  
  cromwellAuthenticationScheme = "service_account"

  // If cromwellAuthenticationScheme is "user_account"
  userAuth {
    // user = ""
    // secretsFile = ""
    // dataStoreDir = ""
  }

  // If cromwellAuthenticationScheme is "service_account"
  serviceAuth {
    pemFile = "/path/to/secret/cromwell-svc-acct.pem"
    serviceAccountId = "806222273987-gffklo3qfd1gedvlgr55i84cocjh8efa@developer.gserviceaccount.com"
  }
}
```

### Data Localization

Data localization can be performed on behalf of an other entity (typically a user).

This allows cromwell to localize file that otherwise wouldn't be accessible using whichever `cromwellAuthenticationScheme` has been defined in the `google` configuration (e.g. if data has restrictive ACLs).
To enable this feature, two pieces of configuration are needed:

**1 - ClientID/Secret**

An entry must be added in the `google` stanza, indicating a pair of client ID / client Secret that have been used to generate a refresh token for the entity that will be used during localization:

```hocon
google {
  cromwellAuthenticationScheme = "service_account"

  serviceAuth {
    pemFile = "/path/to/secret/cromwell-svc-acct.pem"
    serviceAccountId = "806222273987-gffklo3qfd1gedvlgr55i84cocjh8efa@developer.gserviceaccount.com"
  }
  
  userAuthenticationScheme = "refresh"
  
  refreshTokenAuth = {
    client_id = "myclientid.apps.googleusercontent.com"
    client_secret = "clientsecretpassphrase"
  }
}
```

**2 - Refresh Token**

A **refresh_token** field must be specified in the [workflow options](#workflow-options) when submitting the job.  Omitting this field will cause the workflow to fail.

The refresh token is passed to JES along with the client ID and Secret pair, which allows JES to localize and delocalize data as the entity represented by the refresh token.
Note that upon generation of the refresh token, the application must ask for GCS read/write permission using the appropriate scope.

### Docker

It is possible to reference private docker images in dockerhub to be run on JES.
However, in order for the image to be pulled, the docker credentials with access to this image must be provided in the configuration file.

```
docker {
  dockerAccount = "mydockeraccount@mail.com"
  dockerToken = "mydockertoken"
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

The output of this script will be written to a `monitoring.log` file that will be available in the call gcs bucket when the call completes.

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
    defaultZones: "US_Metro MX_Metro"
    defaultDisks: "Disk1 3 SSD, Disk2 500 HDD"
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
| cpu                  |       |   x   |       |
| defaultDisks         |       |   x   |       |
| defaultZones         |       |   x   |       |
| docker               |   x   |   x   |   x   |
| failOnStderr         |   x   |   x   |   x   |
| memory               |       |   x   |       |
| preemptible          |       |   x   |       |

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

```
runtime {
  cpu: 2
}
```

Defaults to "1".

## defaultDisks

Passed to JES: "Disks to attach."

The disks are specified as a comma separated list of disks. Each disk is further separated as a space separated triplet of:

1. Disk name
2. Disk size in GB (not applicable for LOCAL)
3. Disk type

The Disk type must be one of "LOCAL", "SSD", or "HDD". When set to "LOCAL", the size of the drive is automatically provisioned by Google. All disks are set to auto-delete after the job completes.

... more to come ...

```
runtime {
  defaultDisks: "local-disk LOCAL, Disk1 3 SSD, Disk2 500 HDD"
}
```

To change the size of the local disk, set the type of the disk named "local-disk" to a persistent type, and specify the size in GB.

```
runtime {
  defaultDisks: "local-disk 11 SSD"
}
```

Defaults to "local-disk LOCAL".

## defaultZones

Passed to JES: "List of Google Compute Engine availability zones to which resource creation will restricted."

The zones are specified as a space separated list, with no commas.

```
runtime {
  defaultZones: "US_Metro MX_Metro"
}
```

Defaults to "us-central1-a".

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

The memory size is specified as an amount and units of memory, for example "4 G".

```
runtime {
  memory: "4G"
}
```

Defaults to "2G".

## preemptible

Passed to JES: "If applicable, preemptible machines may be used for the run."

... more to come ...

```
runtime {
  preemptible: true
}
```

Defaults to "false".

# Logging

Cromwell accepts three Java Properties for controlling logging:

* `LOG_MODE` - Accepts either `pretty` or `standard` (default `pretty`).  In `standard` mode, logs will be written without ANSI escape code coloring, with a layout more appropriate for server logs, versus `pretty` that is easier to read for a single workflow run.
* `LOG_LEVEL` - Level at which to log (default `info`).
* `LOG_ROOT` - Specifies the directory where logs will be written (default `.`). Currently unused, as logs are only written to standard out, but will be restored in a future update.

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

* **write_to_cache** - Accepts values `true` or `false`.  If `false`, the completed calls from this workflow will not be added to the cache.  See the [Call Caching](#call-caching) section for more details.
* **read_from_cache** - Accepts values `true` or `false`.  If `false`, Cromwell will not search the cache when invoking a call (i.e. every call will be executed unconditionally).  See the [Call Caching](#call-caching) section for more details.
* **outputs_path** - Specifies a path where final workflow outputs will be written.  If this is not specified, workflow outputs will not be copied out of the Cromwell workflow execution directory/path.
* **jes_gcs_root** - (JES backend only) Specifies where outputs of the workflow will be written.  Expects this to be a GCS URL (e.g. `gs://my-bucket/workflows`).  If this is not set, this defaults to the value within `backend.jes.baseExecutionBucket` in the [configuration](#configuring-cromwell).
* **google_project** - (JES backend only) Specifies which google project to execute this workflow.
* **refresh_token** - (JES backend only) Only used if `localizeWithRefreshToken` is specified in the [configuration file](#configuring-cromwell).  See the [Data Localization](#data-localization) section below for more details.
* **auth_bucket** - (JES backend only) defaults to the the value in **jes_gcs_root**.  This should represent a GCS URL that only Cromwell can write to.  The Cromwell account is determined by the `google.authScheme` (and the corresponding `google.userAuth` and `google.serviceAuth`)
* **monitoring_script** - (JES backend only) Specifies a GCS URL to a script that will be invoked prior to the WDL command being run.  For example, if the value for monitoring_script is "gs://bucket/script.sh", it will be invoked as `./script.sh > monitoring.log &`.  The value `monitoring.log` file will be automatically de-localized.

# Call Caching

Call Caching allows Cromwell to detect when a job has been run in the past so it doesn't have to re-compute results.  Cromwell searches the cache of previously run jobs for a one that has the exact same command and exact same inputs.  If a previously run job is found in the cache, Cromwell will **copy the results** of the previous job instead of re-running it.

Cromwell's call cache is maintained in its database.  For best mileage with call caching, configure Cromwell to [point to a MySQL database](#database) instead of the default in-memory database.  This way any invocation of Cromwell (either with `run` or `server` subcommands) will be able to utilize results from all calls that are in that database.

**Call Caching is disabled by default.**  Once enabled, Cromwell will search the call cache for every `call` statement invocation, assuming `read-from-cache` is enabled (see below):

* If there was no cache hit, the `call` will be executed as normal.  Once finished it will add itself to the cache, assuming `read-from-cache` is enabled (see below)
* If there was a cache hit, outputs are **copied** from the cached job to the new job's output directory

> **Note:** If call caching is enabled, be careful not to change the contents of the output directory for any previously run job.  Doing so might cause cache hits in Cromwell to copy over modified data and Cromwell currently does not check that the contents of the output directory changed.

To enable Call Caching, add the following to your Cromwell [configuration](#configuring-cromwell):

```
call-caching {
  enabled = true
  lookup-docker-hash = false
}
```

When `call-caching.enabled=true`, Cromwell will add completed calls to the cache as well as do a cache lookup before running any call.

When `call-caching.lookup-docker-hash=true`, Cromwell will contact external services like DockerHub or Google Container Registry to resolve Docker floating container identifiers like `ubuntu:latest` into immutable hashes while computing the hash of the call invocation.  If this option is false, then the raw value specified in the WDL file for the Docker image is the value that will be used.

Cromwell also accepts two [workflow option](#workflow-options) related to call caching:

* If call caching is enabled, but one wishes to run a workflow but not add any of the calls into the call cache when they finish, the `write_to_cache` option can be set to `false`.  This value defaults to `true`.
* If call caching is enabled, but you don't want to check the cache for any `call` invocations, set the option `read_from_cache` to `false`.  This value also defaults to `true`

> **Note:** If call caching is disabled, the to workflow options `read_from_cache` and `write_to_cache` will be ignored and the options will be treated as though they were 'false'.

# REST API

The `server` subcommand on the executable JAR will start an HTTP server which can accept WDL files to run as well as check status and output of existing workflows.

The following sub-sections define which HTTP Requests the web server can accept and what they will return.  Example HTTP requests are given in [HTTPie](https://github.com/jakubroztocil/httpie) and [cURL](http://curl.haxx.se/)

## REST API Versions

All web server requests include an API version in the url. The current version is `v1`.

## POST /api/workflows/:version

This endpoint accepts a POST request with a `multipart/form-data` encoded body.  The form fields that may be included are:

* `wdlSource` - *Required* Contains the WDL file to submit for execution.
* `workflowInputs` - *Optional* JSON file containing the inputs.  A skeleton file can be generated from the CLI with the [inputs](#inputs) sub-command.
* `workflowOptions` - *Optional* JSON file containing options for this workflow execution.  See the [run](#run) CLI sub-command for some more information about this.

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

workflow sfrazer {
  call x
  call y
  call z {
    input: example="gs://my-bucket/cromwell-executions/sfrazer/example.txt", int=3000
  }
}

--f3fd038395644de596c460257626edd7
Content-Disposition: form-data; name="workflowInputs"; filename="jes0.json"

{
  "sfrazer.x.x": "100"
}

--f3fd038395644de596c460257626edd7
Content-Disposition: form-data; name="workflowOptions"; filename="options.json"

{
  "jes_gcs_root": "gs://sfrazer-dev/workflows"
}

--f3fd038395644de596c460257626edd7--
```

## GET /api/workflows/:version/query

This endpoint allows for querying workflows based on the following criteria:

* `name`
* `status`
* `start` (start datetime)
* `end` (end datetime)

Names and statuses can be given multiple times to include workflows with any of the specified names or statuses.
Valid statuses are `Submitted`, `Running`, `Aborting`, `Aborted`, `Failed`, and `Succeeded`.  `start` and `end` should
be in [ISO8601 datetime](https://en.wikipedia.org/wiki/ISO_8601) format and `start` cannot be after `end`.

cURL:

```
$ curl "http://localhost:8000/api/workflows/v1/query?start=2015-11-01&end=2015-11-03&status=Failed&status=Succeeded"
```

HTTPie:

```
$ http "http://localhost:8000/api/workflows/v1/query?start=2015-11-01&end=2015-11-03&status=Failed&status=Succeeded"
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

## GET /api/workflows/:version/:id/outputs/:call

cURL:

```
$ curl http://localhost:8000/api/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/outputs/three_step.wc
```

HTTPie:

```
$ http http://localhost:8000/api/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/outputs/three_step.wc
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
## GET /api/workflows/:version/:id/logs/:call

This will return paths to the standard out and standard error files that were generated during the execution of a particular fully-qualified name for a call.

A call has one or more standard out and standard error logs, depending on if the call was scattered or not. In the latter case, one log is provided for each instance of the call that has been run.

cURL:

```
$ curl http://localhost:8000/api/workflows/v1/b3e45584-9450-4e73-9523-fc3ccf749848/logs/three_step.wc
```

HTTPie:

```
$ http http://localhost:8000/api/workflows/v1/b3e45584-9450-4e73-9523-fc3ccf749848/logs/three_step.wc
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
        "three_step.wc": [
            {
                "stderr": "/home/user/test/b3e45584-9450-4e73-9523-fc3ccf749848/three_step.wc/stderr6126967977036995110.tmp",
                "stdout": "/home/user/test/b3e45584-9450-4e73-9523-fc3ccf749848/three_step.wc/stdout6128485235785447571.tmp"
            }
        ]
    }
}
```

In the case that the call is inside a `scatter` block, the output for this API will contain a list of stdout/stderr files, one for each shard.  Consider this example:

```
task add_one {
  Int n
  command {
    python -c "print(${n}+1)"
  }
  output {
    Int incremented = read_int(stdout())
  }
}

workflow test {
  Array[Int] list = [1,2,3,4]
  scatter (x in list) {
    call add_one {input: n=x}
  }
}
```

Running this workflow then issuing this API call would return:

```
HTTP/1.1 200 OK
Content-Length: 1256
Content-Type: application/json; charset=UTF-8
Date: Fri, 04 Sep 2015 12:22:45 GMT
Server: spray-can/1.3.3

{
    "id": "cbdefb0f-29ae-475b-a42c-90403f8ff9f8",
    "logs": {
        "test.add_one": [
            {
                "stderr": "/home/user/test/cbdefb0f-29ae-475b-a42c-90403f8ff9f8/call-add_one/shard-0/stderr",
                "stdout": "/home/user/test/cbdefb0f-29ae-475b-a42c-90403f8ff9f8/call-add_one/shard-0/stdout"
            },
            {
                "stderr": "/home/user/test/cbdefb0f-29ae-475b-a42c-90403f8ff9f8/call-add_one/shard-1/stderr",
                "stdout": "/home/user/test/cbdefb0f-29ae-475b-a42c-90403f8ff9f8/call-add_one/shard-1/stdout"
            },
            {
                "stderr": "/home/user/test/cbdefb0f-29ae-475b-a42c-90403f8ff9f8/call-add_one/shard-2/stderr",
                "stdout": "/home/user/test/cbdefb0f-29ae-475b-a42c-90403f8ff9f8/call-add_one/shard-2/stdout"
            },
            {
                "stderr": "/home/user/test/cbdefb0f-29ae-475b-a42c-90403f8ff9f8/call-add_one/shard-3/stderr",
                "stdout": "/home/user/test/cbdefb0f-29ae-475b-a42c-90403f8ff9f8/call-add_one/shard-3/stdout"
            }
        ]
    }
}
```

## GET /api/workflows/:version/:id/logs

This returns a similar format as the `/api/workflows/:version/:id/logs/:call` endpoint, except that it includes the logs for ALL calls in a workflow and not just one specific call.

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
Workflow metadata includes submission, start, and end datetimes, as well as status, inputs and outputs.
Call-level metadata includes inputs, outputs, start and end datetime, backend-specific job id,
return code, stdout and stderr.  Date formats are ISO with milliseconds.


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
Content-Length: 6192
{
  "calls": {
    "sc_test.do_prepare": [{
      "executionStatus": "Done",
      "stdout": "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_prepare/stdout",
      "outputs": {
        "split_files": ["/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_prepare/temp_aa", "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_prepare/temp_ab", "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_prepare/temp_ac", "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_prepare/temp_ad"]
      },
      "inputs": {
        "input_file": "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/workflow-inputs/e46345ba-11.txt"
      },
      "returnCode": 0,
      "backend": "Local",
      "end": "2015-10-01T18:17:56.651-04:00",
      "stderr": "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_prepare/stderr",
      "start": "2015-10-01T18:17:56.204-04:00"
    }],
    "sc_test.do_scatter": [{
      "executionStatus": "Done",
      "stdout": "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_scatter/shard-0/stdout",
      "outputs": {
        "count_file": "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_scatter/shard-0/output.txt"
      },
      "inputs": {
        "input_file": "f"
      },
      "returnCode": 0,
      "backend": "Local",
      "end": "2015-10-01T18:17:56.777-04:00",
      "stderr": "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_scatter/shard-0/stderr",
      "start": "2015-10-01T18:17:56.676-04:00"
    }, {
      "executionStatus": "Done",
      "stdout": "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_scatter/shard-1/stdout",
      "outputs": {
        "count_file": "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_scatter/shard-1/output.txt"
      },
      "inputs": {
        "input_file": "f"
      },
      "returnCode": 0,
      "backend": "Local",
      "end": "2015-10-01T18:17:56.819-04:00",
      "stderr": "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_scatter/shard-1/stderr",
      "start": "2015-10-01T18:17:56.676-04:00"
    }, {
      "executionStatus": "Done",
      "stdout": "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_scatter/shard-2/stdout",
      "outputs": {
        "count_file": "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_scatter/shard-2/output.txt"
      },
      "inputs": {
        "input_file": "f"
      },
      "returnCode": 0,
      "backend": "Local",
      "end": "2015-10-01T18:17:56.806-04:00",
      "stderr": "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_scatter/shard-2/stderr",
      "start": "2015-10-01T18:17:56.676-04:00"
    }, {
      "executionStatus": "Done",
      "stdout": "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_scatter/shard-3/stdout",
      "outputs": {
        "count_file": "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_scatter/shard-3/output.txt"
      },
      "inputs": {
        "input_file": "f"
      },
      "returnCode": 0,
      "backend": "Local",
      "end": "2015-10-01T18:17:56.793-04:00",
      "stderr": "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_scatter/shard-3/stderr",
      "start": "2015-10-01T18:17:56.676-04:00"
    }],
    "sc_test.do_gather": [{
      "executionStatus": "Done",
      "stdout": "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_gather/stdout",
      "outputs": {
        "sum": 11
      },
      "inputs": {
        "input_files": "do_scatter.count_file"
      },
      "returnCode": 0,
      "backend": "Local",
      "end": "2015-10-01T18:17:56.895-04:00",
      "stderr": "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_gather/stderr",
      "start": "2015-10-01T18:17:56.206-04:00"
    }]
  },
  "outputs": {
    "sc_test.do_gather.sum": 11,
    "sc_test.do_prepare.split_files": ["/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_prepare/temp_aa", "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_prepare/temp_ab", "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_prepare/temp_ac", "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_prepare/temp_ad"],
    "sc_test.do_scatter.count_file": ["/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_scatter/shard-0/output.txt", "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_scatter/shard-1/output.txt", "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_scatter/shard-2/output.txt", "/home/jdoe/cromwell/cromwell-executions/sc_test/167dafe7-f474-482a-8194-01a1b770cdfd/call-do_scatter/shard-3/output.txt"]
  },
  "id": "167dafe7-f474-482a-8194-01a1b770cdfd",
  "inputs": {
    "sc_test.do_prepare.input_file": "/home/jdoe/cromwell/11.txt"
  },
  "submission": "2015-10-01T18:17:56.113-04:00",
  "status": "Succeeded",
  "start": "2015-10-01T18:17:56.113-04:00"
}
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

## POST /api/workflows/:version/:id/call-caching

This endpoint allows for reconfiguration of call cache result reuse settings for all calls within a workflow.

Accepted parameters are:

* `allow` Mandatory boolean value, specifies whether call cache result reuse is allowed for all calls in the
   specified workflow.

cURL:

```
$ curl -X POST http://localhost:8000/api/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/call-caching?allow=false
```

HTTPie:

```
$ http POST http://localhost:8000/api/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/call-caching?allow=false
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 17
Content-Type: application/json; charset=UTF-8
Date: Thu, 04 Jun 2015 12:15:33 GMT
Server: spray-can/1.3.3

{
    "updateCount": 3
}

```

## POST /api/workflows/:version/:id/call-caching/:call

This endpoint allows for reconfiguration of call cache result reuse settings for a single call within a workflow.

Accepted parameters are:

* `allow` Mandatory boolean value, specifies whether call cache result reuse is allowed for the specified call in the
  specified workflow.

For scattered calls, individual calls within the scatter can be targeted by appending a dot and the zero-based shard index.
e.g. `scatter_workflow.A.0` would target the zeroth shard of a scattered `A` call.  If a shard index is not supplied for
a scattered call, all shards are targeted for update.

cURL:

```
$ curl -X POST http://localhost:8000/api/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/call-caching/three_step.wc?allow=false
```

HTTPie:

```
$ http POST http://localhost:8000/api/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/call-caching/three_step.wc?allow=false
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 17
Content-Type: application/json; charset=UTF-8
Date: Thu, 04 Jun 2015 12:15:33 GMT
Server: spray-can/1.3.3

{
    "updateCount": 1
}

```

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
