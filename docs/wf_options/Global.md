The following workflow options apply to all workflows and their calls, regardless of the backend being used.


|Option|Values|Description|Example|
|---|---|---|---|
|`write_to_cache`|`true` or `false`|If `false`, the completed calls from this workflow will not be added to the cache.  See the [Call Caching](CallCaching) section for more details.|`write_to_cache: true`|


- **write_to_cache** - Accepts values `true` or `false`.  If `false`, the completed calls from this workflow will not be added to the cache.  See the [Call Caching](CallCaching) section for more details.
    * **read_from_cache** - Accepts values `true` or `false`.  If `false`, Cromwell will not search the cache when invoking a call (i.e. every call will be executed unconditionally).  See the [Call Caching](CallCaching) section for more details.
- **final_workflow_log_dir** - Specifies a path where per-workflow logs will be written.  If this is not specified, per-workflow logs will not be copied out of the Cromwell workflow log temporary directory/path before they are deleted.
- **final_workflow_outputs_dir** - Specifies a path where final workflow outputs will be written.  If this is not specified, workflow outputs will not be copied out of the Cromwell workflow execution directory/path.
- **final_call_logs_dir** - Specifies a path where final call logs will be written.  If this is not specified, call logs will not be copied out of the Cromwell workflow execution directory/path.
- **default_runtime_attributes** - A JSON object where the keys are [Runtime Attributes](RuntimeAttributes) and the values are defaults that will be used through the workflow invocation.  Individual tasks can choose to override these values.  See the [Runtime Attributes](RuntimeAttributes) section for more information.
- **continueOnReturnCode** - Can accept a boolean value or a comma separated list of integers in a string.  Defaults to false.  If false, then only return code of 0 will be acceptable for a task invocation.  If true, then any return code is valid.  If the value is a list of comma-separated integers in a string, this is interpreted as the acceptable return codes for this task.
- **workflow_failure_mode** - What happens after a task fails. Choose from:
        * **ContinueWhilePossible** - continues to start and process calls in the workflow, as long as they did not depend on the failing call
        * **NoNewCalls** - no *new* calls are started but existing calls are allowed to finish
        * The default is `NoNewCalls` but this can be changed using the `workflow-options.workflow-failure-mode` configuration option.
- **backend** - Override the default backend specified in the Cromwell configuration for this workflow only.
