The configuration files are in
[Hocon](https://github.com/typesafehub/config/blob/master/HOCON.md#hocon-human-optimized-config-object-notation).

To create your own configuration file, create a new text file and add your custom configuration. At the start of the
file, include the file `application.conf` at the top before your custom configurations.

```hocon
# include the application.conf at the top
include required(classpath("application"))
```

From there, add other configuration values and/or stanzas with your customizations.

```hocon
# include the application.conf at the top
include required(classpath("application"))

# Add customizations
webservice.port = 58000
```

Your configuration file can specify configuration as JSON-like stanzas like:

```hocon
include required(classpath("application"))
webservice {
  port = 8000
  interface = 0.0.0.0
  binding-timeout = 5s
  instance.name = "reference"
}
```

Or, alternatively, as dot-separated values:

```hocon
include required(classpath("application"))
webservice.port = 8000
webservice.interface = 0.0.0.0
webservice.binding-timeout = 5s
webservice.instance.name = "reference"
```

This allows any value to be overridden on the command line:

```
java -Dwebservice.port=8080 cromwell.jar ...
```

To customize configuration it is recommended that one copies relevant stanzas from `cromwell.examples.conf` into a new
file, modify it as appropriate, then pass it to Cromwell via:

```
java -Dconfig.file=/path/to/yourOverrides.conf cromwell.jar ...
```

A description of options and example stanzas may be found in the file
[`cromwell.examples.conf`](cromwell.examples.conf).

## I/O

Cromwell centralizes as many of its I/O operations as possible through a unique entry point. This allows users to effectively control and throttle the number of requests and resources allocated to those operations throughout the entire system.
It is possible to configure this throttling behavior in the configuration:

```
system.io {
  number-of-requests = 100000
  per = 100 seconds
}
```

This is particularly useful when running Cromwell on a JES backend for example, as Google imposes a quota on the number of GCS queries that can be made.

### Resilience

I/O operations can fail for a number of reason from network failures to server errors. Some of those errors are not fatal and can be retried.
Cromwell will retry I/O operations on such retryable errors, up to a number of times. This number (more precisely the number of attempts that will be made) can be set using the following configuration option:

```
system.io {
  # Number of times an I/O operation should be attempted before giving up and failing it.
  number-of-attempts = 5
}
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

By default batch inserts will be processed in blocks of 2000. To modify this value add the field `insert-batch-size` to the `database` stanza.

Cromwell stores metadata about each job and workflow intended for end users. This metadata includes paths to job
results, start and end times, etc. The metadata grows at a significantly faster rate than the rest of the internal
engine data. To use a separate database for metadata, under the `database` config section, configure a sub-path for
`metadata` with custom settings.

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

If no override is found for `metadata`, Cromwell falls back to using the settings under the root `database`
configuration. This feature should be considered experimental and likely to change in the future.

## SIGINT abort handler

For backends that support aborting task invocations, Cromwell can be configured to automatically try to abort all currently running calls (and set their status to `Aborted`) when a SIGINT is sent to the Cromwell process.  To turn this feature on, set the configuration option

```
system {
  abort-jobs-on-terminate=true
}
```

Or, via `-Dsystem.abort-jobs-on-terminate=true` command line option.

By default, this value is false when running `java -jar cromwell.jar server`, and true when running `java -jar cromwell.jar run <workflow source> <inputs>`.