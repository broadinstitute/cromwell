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

... to be determined ...
