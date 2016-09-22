# Cromwell Change Log

## 0.20

* The default per-upload bytes size for GCS is now the minumum 256K
instead of 64M. There is also an undocumented config key
`google.upload-buffer-bytes` that allows adjusting this internal value.

* Updated Docker Hub hash retriever to parse json with [custom media
types](https://github.com/docker/distribution/blob/05b0ab0/docs/spec/manifest-v2-1.md).

* Added a `/batch` submit endpoint that accepts a single wdl with
multiple input files.

* The `/query` endpoint now supports querying by `id`, and submitting
parameters as a HTTP POST.

## 0.21


* Warning: Significant database updates when you switch from version 0.19 to 0.21 of Cromwell.
  Please be refer to MIGRATION.md for more details.

* There are significant architectural changes related related to increases in performance and scaling.

* The biggest changes from 0.19 to 0.21 are related to the application.conf file, which has been restructured significantly.
The configuration for backends now is all contained within a `backend` stanza, which specifies 1 stanza per name per backend and a default backend, as follows:

```
backend {
    default=Local
    providers {
        Local {
            actor-factory: "class path to BackendLifecycleActorFactory implementation"
            config {
                ... backend specific config ...
            }
        }
        JES {
            actor-factory: "class path to BackendLifecycleActorFactory implementation"
            config {
                ... backend specific config ...
            }
        }
        SGE {
            actor-factory: "class path to BackendLifecycleActorFactory implementation"
            config {
                ... backend specific config ...
            }
        }
    }
}
```
* A new `/stats` endpoint has been added to get workflow and job count for a Cromwell running in server mode.

* Renamed Workflow Options:
   “workflow_log_dir” -> “final_workflow_log_dir”
    “call_logs_dir” -> “final_call_logs_dir”
    “outputs_path” -> “final_workflow_outputs_dir”
    “defaultRuntimeOptions” -> “default_runtime_attributes”

* Timing diagrams endpoint has been updated to include additional state information about jobs.

* Add support for Google Private IPs through `noAddress` runtime attribute. If set to true, the VM will NOT be provided with a public IP address.
*Important*: Your project must be whitelisted in "Google Access for Private IPs Early Access Program". If it's not whitelisted and you set this attribute to true, the task will hang.
  Defaults to `false`.
  e.g:
```
task {
    command {
        echo "I'm private !"
    }
    
    runtime {
        docker: "ubuntu:latest"
        noAddress: true
    }
}
```

* The Local and the SGE backend have been merged into a generic
Shared File System (SFS) backend. This updated backend can be configured
to work with various other command line dispatchers such as LSF. See the
[README](README.md#sun-gridengine-backend) for more info.

* On the JES and SFS backends, task `command` blocks are now always
passed absolute paths for input `File`s.

* On the SFS backends, the call directory now contains two sub-directories:
    * `inputs` contains all the input files that have been localized for this task (see next below for more details)
    * `execution` contains all other files (script, logs, rc, potential outputs etc...)

* Override the default database configuration by setting the keys
`database.driver`, `database.db.driver`, `database.db.url`, etc.
* Override the default database configuration by setting the keys
`database.driver`, `database.db.driver`, `database.db.url`, etc. 

For example:
```
# use a mysql database
database {
  driver = "slick.driver.MySQLDriver$"
  db {
    driver = "com.mysql.jdbc.Driver"
    url = "jdbc:mysql://host/cromwell"
    user = "user"
    password = "pass"
    connectionTimeout = 5000
  }
}
```

