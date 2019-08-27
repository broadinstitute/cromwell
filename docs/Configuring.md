## Overview

You can configure Cromwell settings either through configuration files or the Java command line.

Check out the tutorial on [How to Configure Cromwell](tutorials/ConfigurationFiles) for more information.

### Configuration examples

You can find a description of options and example stanzas in the [Cromwell Example Configuration][cromwell-examples-conf],
along with backend provider examples in the [Example Providers Folder][cromwell-examples-folder].

### Custom configuration files

You write configuration files in
[HOCON](https://github.com/typesafehub/config/blob/master/HOCON.md#hocon-human-optimized-config-object-notation).

To run using your configuration file, you should copy relevant stanzas from `cromwell.examples.conf` into a new
file, modify it as appropriate, then pass it to Cromwell via:

```
$ java -Dconfig.file=/path/to/yourOverrides.conf cromwell.jar ...
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

_JSON-like stanza:_

```hocon
include required(classpath("application"))
webservice {
  port = 8000
  interface = 0.0.0.0
}
```

_Dot-separated values:_

```hocon
include required(classpath("application"))
webservice.port = 8000
webservice.interface = 0.0.0.0
```

## Configuration via command line

In addition to using configuration files, you can use dot-separated configuration names to specify values directly on the Java command line:

```
$ java -Dwebservice.port=8080 cromwell.jar ...
```

## Advanced configuration

**WARNING:** These advanced configuration values can significantly affect the performance of Cromwell. 

### Server

By default the Cromwell server will bind to `0.0.0.0` on port `8000`.  
You can then access it through a browser at `http://localhost:8000`.  
To change these settings, simply edit the following values in your configuration file:

```
webservice {
  port = 9000
  interface = 0.0.0.0
}
```

The above configuration will use port `9000`.

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

***Abort configuration***

Cromwell will scan for abort requests using default configuration values equivalent to those below. In most circumstances
there shouldn't be a need to override these defaults.

```hocon
system {
  abort {
    # How frequently Cromwell should scan for aborts.
    scan-frequency: 30 seconds

    # The cache of in-progress aborts. Cromwell will add entries to this cache once a WorkflowActor has been messaged to abort.
    # If on the next scan an 'Aborting' status is found for a workflow that has an entry in this cache, Cromwell will not ask
    # the associated WorkflowActor to abort again.
    cache {
      # Guava cache concurrency.
      concurrency: 1
      # How long entries in the cache should live from the time they are added to the cache.
      ttl: 20 minutes
      # Maximum number of entries in the cache.
      size: 100000
    }
  }
}
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
    driver = "com.mysql.cj.jdbc.Driver"
    url = "jdbc:mysql://host/cromwell?rewriteBatchedStatements=true"
    user = "user"
    password = "pass"
    connectionTimeout = 5000
  }
}
```

To see the full list of possible parameters and values for the `db` stanza see [the slick documentation](http://slick.lightbend.com/doc/3.2.0/api/index.html#slick.jdbc.JdbcBackend$DatabaseFactoryDef@forConfig(String,Config,Driver):Database).

**Cromwell server on MySQL Database**

You can use [docker-compose](https://github.com/broadinstitute/cromwell/tree/develop/scripts) to link together a Cromwell docker image (built locally with `sbt docker` or available on [Dockerhub](https://hub.docker.com/r/broadinstitute/cromwell/)) with a MySQL docker image.

To change the version of Cromwell used, [change the tag in `compose/cromwell/Dockerfile`](https://github.com/broadinstitute/cromwell/blob/develop/scripts/docker-compose-mysql/compose/cromwell/Dockerfile).

**Local**

`docker-compose up` from this directory will start a Cromwell server running on a MySQL instance with local backend.

The default configuration file used can be [found at `compose/cromwell/app-config/application.conf`](https://github.com/broadinstitute/cromwell/blob/develop/scripts/docker-compose-mysql/compose/cromwell/app-config/application.conf).
To override it, simply mount a volume containing your custom `application.conf` to `/app-config` ([see `jes-cromwell/docker-compose.yml` for an example](https://github.com/broadinstitute/cromwell/blob/develop/scripts/docker-compose-mysql/jes-cromwell/docker-compose.yml)).

**Google Cloud**

The [`jes-cromwell` directory](https://github.com/broadinstitute/cromwell/tree/develop/scripts/docker-compose-mysql/jes-cromwell) is an example of how to customize the original compose file with a configuration file and environment variables.

It uses the application default credentials of the host machine. To use it make sure your gcloud is up to date and that your [application-default credentials](https://developers.google.com/identity/protocols/application-default-credentials) are set up.

Then run `docker-compose -f docker-compose.yml -f jes-cromwell/docker-compose.yml up` to start a Cromwell server with a Google Cloud backend on MySQL.

**MySQL**

The data directory in the MySQL container is [mounted to `compose/mysql/data`](https://github.com/broadinstitute/cromwell/tree/develop/scripts/docker-compose-mysql/compose/mysql/init), which allows the data to survive a `docker-compose down`.

To disable this feature, simply remove the `./compose/mysql/data:/var/lib/mysql` line in the [volume section of `docker-compose.yml`](https://github.com/broadinstitute/cromwell/blob/develop/scripts/docker-compose-mysql/docker-compose.yml).

Note that in such case, the data will still be preserved by a `docker-compose stop` that stops the container but doesn't delete it.

**Notes**

To run Cromwell in the background, add `-d` at the end of the command:
`docker-compose up -d`.

To then see the logs for a specific service, run `docker-compose logs -f <service>`. 
For example `docker-compose logs -f cromwell`.

For more information about docker compose: [Docker compose doc](https://docs.docker.com/compose/).

**Insert Batch Size**

Cromwell queues up and then inserts batches of records into the database for increased performance. You can adjust the
number of database rows batch inserted by Cromwell as follows:

```hocon
database {
  insert-batch-size = 2000
}
```

**Separate Metadata Database**

This feature should be considered _experimental_ and likely to change in the future.

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

**Database Time Zones**

Cromwell's default configuration assumes that its MySQL database is set to UTC.

The following MySQL configurations typically default to UTC and work with Cromwell out of the box:
- Google CloudSQL
- An official MySQL image running in Docker

These configurations may use the system, or local, time zone instead:
- MySQL installed natively on a workstation or server

If Cromwell fails to start with a message like
```
The server time zone value 'XXX' is unrecognized or represents more than one time zone.
```
you can resolve the problem by adding the option `&serverTimezone=UTC` to your database connection URL:
```hocon
url = "jdbc:mysql://host/cromwell?rewriteBatchedStatements=true&serverTimezone=UTC"
```

Using this option does not alter your database's underlying timezone; rather, it causes Cromwell to "speak UTC" when communicating with the DB, and the DB server performs the conversion for you. 

**Using Cromwell with Postgresql**

To use Postgresql as the database, you will need to install and enable the
Large Object extension.  If the extension is present, setting up the database
requires just these commands:

```
$ createdb cromwell
$ psql -d cromwell -c "create extension lo;"
```

Postgresql configuration in Cromwell is very similar to MySQL.  An example:

```hocon
database {
  profile = "slick.jdbc.PostgresProfile$"
  db {
    driver = "org.postgresql.Driver"
    url = "jdbc:postgresql://localhost:5432/cromwell"
    user = "user"
    password = "pass"
    port = 5432
    connectionTimeout = 5000
  }
}
```

### Abort

**Control-C (SIGINT) abort handler**

For backends that support aborting jobs, Cromwell can be configured to automatically try to abort all calls when it receives a Control-C, also known as SIGINT. All currently running calls will also set their status to `Aborted`.

To explicitly turn this feature on or off, set the configuration option:

```hocon
system {
  abort-jobs-on-terminate=true
}
```

Or, via `-Dsystem.abort-jobs-on-terminate=true` command line option.

By default, this value is false when running `java -jar cromwell.jar server`, and true when running `java -jar cromwell.jar run <workflow source> <inputs>`.

Read the [Abort](execution/ExecutionTwists/#abort) section to learn more about how abort works.

### Call caching

Call Caching allows Cromwell to detect when a job has been run in the past so it doesn't have to re-compute results.  To learn more see [Call Caching](cromwell_features/CallCaching).

To enable Call Caching, add the following to your Cromwell configuration:

```
call-caching {
  enabled = true
  invalidate-bad-cache-results = true
}
```

When `call-caching.enabled=true` (default: `false`), Cromwell will be able to to reference or copy results from previously run jobs (when appropriate).
When `invalidate-bad-cache-results=true` (default: `true`), Cromwell will invalidate any cache results which contain files that cannot be accessed within a cache-hit. This is usually desired, but might be unwanted if this failure occurs for external reasons, such as a difference in user authentication.

Cromwell also accepts [Workflow Options](wf_options/Overview#call-caching-options) to override the cache read/write behavior.  

### Local filesystem options

When running a job on the Config (Shared Filesystem) backend, Cromwell provides some additional options in the backend's config section:

```HOCON
      config {
        filesystems {
          local {
            caching {
              # When copying a cached result, what type of file duplication should occur. Attempted in the order listed below:
              duplication-strategy: [
                "hard-link", "soft-link", "copy"
              ]

              # Possible values: file, path, path+modtime
              # "file" will compute an md5 hash of the file content.
              # "path" will compute an md5 hash of the file path. This strategy will only be effective if the duplication-strategy (above) is set to "soft-link",
              # in order to allow for the original file path to be hashed.
              # "path+modtime" will compute an md5 hash of the file path and the last modified time. The same conditions as for "path" apply here.
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

### Workflow log directory

To change the directory where Cromwell writes workflow logs, change the directory location via the setting:

```hocon
workflow-options {
    workflow-log-dir = "cromwell-workflow-logs"
}
```


**Preserving Workflow Logs**

By default Cromwell erases the per workflow logs when the workflow completes to reduce disk usage. You can change this behavior by setting the following value to `false`:

```hocon
workflow-options {
    workflow-log-temporary = true
}
```


**Exception monitoring via Sentry**

Cromwell supports [Sentry](https://docs.sentry.io) which is a service that can be used to monitor exceptions reported in an application’s logs.

To enable Sentry monitoring in Cromwell, enter your DSN URL using the system property:

```properties
sentry.dsn=DSN_URL
```


### Job shell configuration

Cromwell allows for system-wide or per-backend job shell configuration for running user commands rather than always
using the default `/bin/bash`. To set the job shell on a system-wide basis use the configuration key `system.job-shell` or on a
per-backend basis with `<config-key-for-backend>.job-shell`. For example:

```
# system-wide setting, all backends get this
-Dsystem.job-shell=/bin/sh
```

```
# override for just the Local backend
-Dbackend.providers.Local.config.job-shell=/bin/sh
```

For the Config backend the value of the job shell will be available in the `${job_shell}` variable. See Cromwell's `reference.conf` for an example
of how this is used for the default configuration of the `Local` backend.

[cromwell-examples-conf]: https://www.github.com/broadinstitute/cromwell/tree/develop/cromwell.examples.conf
[cromwell-examples-folder]: https://www.github.com/broadinstitute/cromwell/tree/develop/cromwell.example.backends

### Workflow Heartbeats

**Cromwell ID**

Each Cromwell instance is given a `cromwell_id` that is either randomly generated or configured.

By default, the Cromwell ID is `cromid-<7_digit_random_hex>`.

A custom identifier may replace the "cromid" portion of the string. For example:

```hocon
system {
  cromwell_id = "main"
}
```

This would generates a `cromwell_id` of `main-<7_digit_random_hex>`. Each time Cromwell restarts the random part of the
ID will change, however the `main` prefix would remain the same.

If the random part of the Cromwell ID should not be generated, set the configuration value:

```hocon
system {
  cromwell_id_random_suffix = false
}
```

**Heartbeat TTL**

When a Cromwell instance begins running or resuming a workflow it stores the above `cromwell_id` within the database row
for the workflow, along with a timestamp called the "heartbeat". As the workflow continues to run the Cromwell instance
will intermittently update the heartbeat for the running workflow.  If the Cromwell dies, after some time-to-live (TTL),
the workflow has been abandoned, and will be resumed by another available Cromwell instance.

Adjust the heartbeat TTL via the configuration value:

```hocon
system.workflow-heartbeats {
  ttl = 10 minutes
}
```

The default TTL is 10 minutes. The shortest the TTL option is 10 seconds.

**Heartbeat Interval**

The interval for writing heartbeats may be adjusted via:

```hocon
system.workflow-heartbeats {
  heartbeat-interval = 2 minutes
}
```

The default interval is 2 minutes. The shortest interval option is 3.333 seconds. The interval may not be greater than
the TTL.

**Heartbeat Failure Shutdown**

Cromwell will automatically shutdown when unable to write heartbeats for a period of time. This period of time may be
adjusted via:

```hocon
system.workflow-heartbeats {
  write-failure-shutdown-duration = 5 minutes
}
```

The default shutdown duration is 5 minutes. The maximum allowed shutdown duration is the TTL.

**Heartbeat Batch Size**

Workflow heartbeats are internally queued by Cromwell and written in batches. When the configurable batch size is
reached, all of the heartbeats within the batch will be written at the same time, even if the heartbeat interval has not
elapsed.

This batch threshold may be adjusted via:

```hocon
system.workflow-heartbeats {
  write-batch-size = 100
}
```

The default batch size is 100.

**Heartbeat Threshold**

Cromwell writes one batch of workflow heartbeats at a time. While the internal queue of heartbeats-to-write passes above
a configurable threshold then [instrumentation](developers/Instrumentation.md) may send a metric signal that the
heartbeat load is above normal.

This threshold may be configured via the configuration value:

```hocon
system.workflow-heartbeats {
  write-threshold = 100
}
```

The default threshold value is 100, just like the default for the heartbeat batch size.

### YAML

**Maximum number of nodes**

Cromwell will throw an error when detecting cyclic loops in Yaml inputs. However one can craft small acyclic YAML
documents that consume significant amounts of memory or cpu. To limit the amount of processing during parsing, there is
a limit on the number of nodes parsed per YAML document.

This limit may be configured via the configuration value:

```hocon
yaml {
  max-nodes = 1000000
}
```

The default limit is 1,000,000 nodes.

**Maximum nesting depth**

There is a limit on the maximum depth of nested YAML. If you decide to increase this value, you will likely need to also
increase the Java Virtual Machine's thread stack size as well using
[either `-Xss` or `-XX:ThreadStackSize`](https://docs.oracle.com/javase/8/docs/technotes/tools/unix/java.html).

This limit may be configured via the configuration value:

```hocon
yaml {
  max-depth = 1000
}
```

The default limit is a maximum nesting depth of 1,000.
