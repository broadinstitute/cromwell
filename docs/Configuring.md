_For the Doc-A-Thon_  

**NOTES TO EDITOR:**
- There are a large number of options in [`cromwell.examples.conf`](https://github.com/broadinstitute/cromwell/blob/develop/cromwell.examples.conf). Only the existing options listed here were updated. We can add more in the future.
- These subsections are not headers, so links to subsections do not direct you there. (Note on note, there are now some subsections. Others are just **bold**.)

**Questions to answer and things to consider:**

1. Who is visiting the Configuration page?  
*Are they looking to make a config file? Edit one?*
2. What do they need to know first?
3. Is all the important information there? If not, add it!  
*What about useful info like example configs, like Google cloud. REMEMBER to take out your personal details!*
4. Are there things that don't need to be there? Remove them.
5. Are the code and instructions accurate? Try it!

---
 **DELETE ABOVE ONCE COMPLETE**

---

## Overview

You can configure cromwell settings either through configuration files or the java command line.

### cromwell.examples.conf

You can find a description of options and example stanzas in the file
[`cromwell.examples.conf`](https://github.com/broadinstitute/cromwell/blob/develop/cromwell.examples.conf).

### Custom Configuration Files

You write configuration files in
[HOCON](https://github.com/typesafehub/config/blob/master/HOCON.md#hocon-human-optimized-config-object-notation).

To run using your configuration file, you should copy relevant stanzas from `cromwell.examples.conf` into a new
file, modify it as appropriate, then pass it to Cromwell via:

```
java -Dconfig.file=/path/to/yourOverrides.conf cromwell.jar ...
``` 

To create your own configuration file, start by creating a new text file, for example `my.conf`.

At the start of your file, include the file `application.conf` at the top before your custom configurations.

```hocon
# include the application.conf at the top
include required(classpath("application"))
```

From there, copy or add other configuration values and/or stanzas with your customizations.

```hocon
# include the application.conf at the top
include required(classpath("application"))

# Add customizations
webservice.port = 58000
```

Your configuration file can specify configuration as JSON-like stanzas or as dot-separated values. These next two examples are are equivalent.

JSON-like stanza:

```hocon
include required(classpath("application"))
webservice {
  port = 8000
  interface = 0.0.0.0
  binding-timeout = 5s
  instance.name = "reference"
}
```

Dot-separated values:

```hocon
include required(classpath("application"))
webservice.port = 8000
webservice.interface = 0.0.0.0
webservice.binding-timeout = 5s
webservice.instance.name = "reference"
```

## Configuration via Command Line Arguments

In addition to using configuration files, you can use dot-separated configuration names to specify values directly on the java command line:

```
java -Dwebservice.port=8080 cromwell.jar ...
```

## Advanced Configuration

**WARNING:** These advanced configuration values can significantly affect the performance of Cromwell.

### Server

By default the Cromwell server will bind to `0.0.0.0` on port `8000`.
It can then be can be accessed through a browser at `http://localhost:8000`.
To change these settings, simply edit the following values in your configuration file:

```
webservice {
  port = 9000
  interface = 0.0.0.0
}
```

The above configuration will for instance use port `9000`.

Cromwell uses `akka-http` to serve requests. For more advanced configuration settings, refer to the [akka-http](https://doc.akka.io/docs/akka-http/current/scala/http/configuration.html) documentation.

For example, to increase the request timeout to 30 seconds you can add this stanza to your configuration file:

```
akka.http.server.request-timeout = 30s
```

### I/O

**I/O Throttling**

Certain [backends](backends/Backends) impose I/O limits. For example the Pipelines API imposes a quota on the number of queries that can be made per second.

You can effectively control and throttle the number of requests and resources allocated to those operations in the `system.io` configuration:

```
system.io {
  number-of-requests = 100000
  per = 100 seconds
}
```

**I/O Resilience**

I/O operations can fail for a number of reason from network failures to server errors. Some of those errors are not fatal and can be retried.

Cromwell will retry I/O operations on such retryable errors, for a limited number of times before giving up and failing. This number (more precisely the number of attempts that will be made) can be set using the following configuration option:

```
system.io {
  number-of-attempts = 5
}
```

### Workflows

**Max Concurrent Workflows**

Cromwell has a configurable cap on the number of workflows running at a time. You can adjust the limit from the default `5000` by setting:

```hocon
system.max-concurrent-workflows = 5000
```

**New Workflow Poll Rate**

Cromwell will look for new workflows to start on a regular interval, configured as a number of seconds. You can change the polling rate from the default `2` seconds by editing the value:

```hocon
system.new-workflow-poll-rate = 2
```

**Max Workflow Launch Count**

On every poll, Cromwell will take at limited number of new submissions, provided there are new workflows to launch and the `system.max-concurrent-workflows` number has not been reached. While the default is to launch up to `50` workflows, you can override this by setting:

```hocon
system.max-workflow-launch-count = 50
```

### Database

**Using a MySQL Database**

Cromwell tracks the execution of workflows and stores outputs of task invocations in a SQL database. Cromwell supports either an external MySQL database, or a temporary in-memory database.

By default, Cromwell uses an in-memory database which will only live for the duration of the JVM.  This provides a quick way to run workflows locally without having to set up MySQL, though it also makes workflow executions somewhat transient.

To configure Cromwell to instead point to a MySQL database, first create the empty database.  In the example below, the database name is `cromwell`.

Then, edit your configuration file `database` stanza, as follows:

```hocon
database {
  profile = "slick.jdbc.MySQLProfile$"
  db {
    driver = "com.mysql.jdbc.Driver"
    url = "jdbc:mysql://host/cromwell?rewriteBatchedStatements=true"
    user = "user"
    password = "pass"
    connectionTimeout = 5000
  }
}
```

**Insert Batch Size**

Cromwell queues up and then inserts batches of records into the database for increased performance. You can adjust the number of database rows batch inserted by cromwell as follows:

```hocon
database {
  insert-batch-size = 2000
}
```

**Separate Metadata Database**

This feature should be considered experimental and likely to change in the future.

Cromwell stores metadata about each job and workflow intended. This metadata is intended for end users, and includes paths to job results, start and end times, etc. The metadata grows at a significantly faster rate than the rest of the internal engine data.

To use a separate database for metadata, under the `database` config section, configure a sub-path for `metadata` with custom settings.

```hocon
database {
  # Store metadata in a file on disk that can grow much larger than RAM limits.
  metadata {
    profile = "slick.jdbc.HsqldbProfile$"
    db {
      driver = "org.hsqldb.jdbcDriver"
      url = "jdbc:hsqldb:file:metadata-db-file-path;shutdown=false;hsqldb.tx=mvcc"
      connectionTimeout = 3000
    }
  }
}
```

If no override is found for `metadata`, Cromwell falls back to using the settings under the root `database` configuration.

## Abort

**Control-C (SIGINT) abort handler**

For backends that support aborting jobs, Cromwell can be configured to automatically try to abort all calls when it receives a Control-C, also known as SIGINT. All currently running calls will also set their status to `Aborted`.

By default, this value is `false` when running `java -jar cromwell.jar server`, and `true` when running `java -jar cromwell.jar run <workflow source> <inputs>`.

To explicitly turn this feature on or off, set the configuration option:

```hocon
system {
  abort-jobs-on-terminate=true
}
```

Or, via `-Dsystem.abort-jobs-on-terminate=true` command line option.

By default, this value is false when running `java -jar cromwell.jar server`, and true when running `java -jar cromwell.jar run <workflow source> <inputs>`.


**Call Caching**

Call Caching allows Cromwell to detect when a job has been run in the past so it doesn't have to re-compute results.  To learn more see [Call Caching](CallCaching).

To enable Call Caching, add the following to your Cromwell configuration:

```
call-caching {
  enabled = true
  invalidate-bad-cache-results = true
}
```

When `call-caching.enabled=true` (default: `false`), Cromwell will be able to to reference or copy results from previously run jobs (when appropriate).
When `invalidate-bad-cache-results=true` (default: `true`), Cromwell will invalidate any cache results which contain files that cannot be accessed within a cache-hit. This is usually desired, but might be unwanted if this failure occurs for external reasons, such as a difference in user authentication.

Cromwell also accepts [Workflow Options](WorkflowOptions) to override the cache read/write behavior.  

**Local Filesystem Options**

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