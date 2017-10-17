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
    * **monitoring_script** - (JES backend only) Specifies a GCS URL to a script that will be invoked prior to the user command being run.  For example, if the value for monitoring_script is "gs://bucket/script.sh", it will be invoked as `./script.sh > monitoring.log &`.  The value `monitoring.log` file will be automatically de-localized.
