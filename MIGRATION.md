0.19 to 0.21
============

## Configuration

The biggest changes from 0.19 to 0.21 are related to the application.conf file, which has been restructured significantly.

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

### Migrating a 0.19 config file to a 0.21 config file

Copy `core/src/main/resources/application.conf` somewhere, and open up your old config file and copy the following values from the 0.19 config into the 0.21 config:

Below is a table based diff, but there is also a [visual diff](http://i.imgur.com/i8j0Zvb.jpg) of the changes.

|0.19|0.21|Comments|
|----|----|--------|
||`akka.actor.deployment`|Configures the number of threads to handle workflow log copying|
||`akka.dispatchers.io-dispatcher`|Configures the "IO" dispatcher|
|`google.applicationName`|`google.application-name`||
|`google.cromwellAuthenticationScheme`|removed|re-expressed as other configuration settings|
|`google.serverAuth.pemFile`|`google.auths` with scheme `service_account`, then `pem-file`||
|`google.serverAuth.serviceAccountId`|`google.auths` with scheme `service_account`, then `service-account-id`||
|`google.userAuthenticationScheme`|re-expressed as other configuration settings||
|`google.refreshTokenAuth.client_id`|`google.auths` with scheme `refresh_token`, then `client-id`||
|`google.refreshTokenAuth.client_secret`|`google.auths` with scheme `refresh_token`, then `client-secret`||
||`engine.filesystems`|Defines additional filesystems that the engine is able to use to execute read and write WDL functions at the workflow-level|
|`docker`|`backend.JES.config.dockerhub`|Moved to JES-specific configuration|
|`jes.project`|`backend.JES.config.project`||
|`jes.baseExecutionBucket`|`backend.JES.config.root`||
|`jes.endpointUrl`|`backend.JES.config.genomics.endpoint-url`||
|`jes.maximumPollingInterval`|`backend.JES.config.maximum-polling-interval`||
|`shared-filesystem.root`|`backend.[Local|SGE].config.root`||
|`shared-filesystem.localization`|`backend.[Local|SGE].config.filesystems.localization`||
|`backend.abortJobsInTerminate`|`system.abort-jobs-on-terminate`||
|`restart-workflow`|`system.workflow-restart`||
|`database.config`|unchanged|Set your database|
|`database.main.[mysql|hsqldb]`|unchanged|Set your database credentials and connection information|

## Workflow Options

The following workflow options have changed names:

|0.19|0.21|
|----|----|
|workflow_log_dir|final_workflow_log_dir|
|outputs_path|final_workflow_outputs_dir|
|call_logs_dir|final_call_logs_dir|
|defaultRuntimeOptions|default_runtime_attributes|
|workflowFailureMode|workflow_failure_mode|

## API Endpoints

* `POST /api/workflows/:version/validate` has been removed
* `GET /api/workflows/:version/:id/outputs/:call` has been removed
* `GET /api/workflows/:version/query` now requires start and end offsets.  See the README section for this API for more details.

## Database Migrations

Cromwell's database structure has changed significantly between 0.19 and 0.21.
All pre-existing data has to be transformed/moved to new tables in order to be usable.
The migration process can be split into 2 steps:

### Restart Migration

Cromwell 0.21 will migrate workflows that were in **Running** state to the new database structure in order to attempt to resume them once the server restarts.
This will ensure that even if cromwell is stopped while workflows are still running they aren't lost.
No particular action/configuration is required for this step.

### Metadata Migration (MySQL Only)

In order to keep metadata from previous (and current) workflow runs, all the data has to be moved to a new centralized table.
Depending on the number and shape of the workflows in your database **this step can be significantly time and space consuming**.
It is not possible to give an accurate estimation due to the multiple variables in play like number of workflows, complexity (number of tasks, scatters, attempts per task, etc...), hardware performance (of the database, of the machine running cromwell), backends used, etc...

However, a good rule of thumb is to make sure that **your database has enough disk space to grow a factor 10**.
This is due to the fact that data is de-normalized during the migration. In particular all inputs and outputs, which means the more complex / large your outputs are (large arrays, etc...), the more your database will grow.
Also be aware that **the migration can take several hours for substantially large databases.**

#### Important Notes

* For better performance, make sure the flag `rewriteBatchedStatements` is set to `true`. This can be done by adding to your database connection url.

e.g:

    jdbc:mysql:http://localhost:3006/cromwell_db?rewriteBatchedStatements=true
    
See the [mysql doc](https://dev.mysql.com/doc/connector-j/5.1/en/connector-j-reference-configuration-properties.html) for more information.

* Because of the nature of wdl inputs and outputs, as well as the way they were stored up until cromwell 0.19, it is necessary to pull them out of the database, process them one by one in cromwell, and re-insert them.
For that reason, if your database contains workflows with very large inputs or outputs (for example `Array`s of several thousands of elements, large matrix, etc...),
you might want to tune the `migration-read-batch-size` and `migration-write-batch-size` configuration fields (see the `database` section in `application.conf`)

    * `migration-read-batch-size` sets the number of rows that should be retrieved at the same time from the database. Once every row is processed, the next set is retrieved until there are no more.
If your workflows/tasks have very large inputs or outputs, a number too large here could cause out of memory errors, or extended waiting times to pull the data from the database.
On the other hand, a number too small could decrease performance by causing more round-trips to the database with the associated overhead.
 
    * `migration-write-batch-size` sets the number of insert statements that are buffered before being committed in a transaction. This decreases the number of queries to the database, while making sure the batch never gets too big.
You might consider decreasing this field if your workflows/tasks have very large `String`s as input or output (for example several MB of text).