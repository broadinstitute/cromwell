# Workflow Options Overview

Workflow options can affect the execution of a single workflow without having to change configuration options or restart Cromwell. 

You provide workflow options to Cromwell in a JSON format. This can be supplied at workflow-submit time either via the [CLI](../CommandLine/) or the [REST endpoint](../api/RESTAPI/):

```json
{
	"option_name_1": "option value 1",
	"option_name_2": "option value 2"
}
```

Unless otherwise specified you can expect workflow options to override any hard-coded defaults in Cromwell or defaults provided in the [configuration file](../Configuring), but to be overridden by any values provided in the workflow definition file itself (WDL or CWL).

Some workflow options apply only to tasks running on the [Google Pipelines API backend](Google).

# Global Workflow Options 

The following workflow options apply to all workflows and their calls regardless of the backend being used.

## Runtime Attributes

Some options allow you to override or set defaults for runtime attributes.

### Setting Default Runtime Attributes

You can supply a default for any [Runtime Attributes](../RuntimeAttributes) by adding a `default_runtime_attributes` map to your workflow options file. Use the key to provide the attribute name and the value to supply the default. 

These defaults replace any defaults in the Cromwell configuration file but are themselves replaced by any values explicitly provided by the task in the WDL or CWL file.

Example `options.json`:
```json
{
    "default_runtime_attributes": {
        "docker": "ubuntu:latest",
        "continueOnReturnCode": [4, 8, 15, 16, 23, 42]
    }
}
```

In this example, if a task in a workflow specifies a `docker:` attribute, the task will get what it specifies. However if any task does not provide a value then it will be treated as though it had specified `ubuntu:latest`.

### Specific Runtime Attributes

|Option|Value|Description|
|---|---|---|
|`continueOnReturnCode`|`true` or `false` or integer array|Globally overrides the `continueOnReturnCode` [runtime attribute](../RuntimeAttributes) for all tasks| 
|`backend`|An [available](../Configuring) backend|Set the **default** backend specified in the Cromwell configuration for this workflow only.|

Example `options.json`:
```json
{
    "continueOnReturnCode": false,
    "backend": "Local"
}
```

In this example, all tasks will be given to the `Local` backend unless they provide a value explicitly in their `runtime { ... }` block. In addition, the `continueOnReturnCode` value for all tasks is hard-coded to `false`, regardless of what the tasks put in their `runtime` block. **TODO or is just a default ala `default_runtime_attributes`?**

## Workflow Failure

The `workflow_failure_mode` option can be given the following values. This overrides any default set by the `workflow-options.workflow-failure-mode` [configuration](../Configuring) options.

|Value|Description|
|---|---|
|`ContinueWhilePossible`|Continues to start and process calls in the workflow, as long as they did not depend on the failing call.|
|`NoNewCalls`|No *new* calls are started but existing calls are allowed to finish.|

Example `options.json`:
```json
{
    "workflow_failure_mode": "ContinueWhilePossible"
}
```

## Output Copying
|Option|Value|Description|
|---|---|---|
|`final_workflow_outputs_dir`|A directory available to Cromwell|Specifies a path where final workflow outputs will be written. If this is not specified, workflow outputs will not be copied out of the Cromwell workflow execution directory/path.|
|`use_relative_output_paths`| A boolean | When set to `true` this will copy all the outputs relative to their execution directory. my_final_workflow_outputs_dir/~~MyWorkflow/af76876d8-6e8768fa/call-MyTask/execution/~~output_of_interest . Cromwell will throw an exception when this leads to collisions. When the option is not set it will default to `false`.|
|`final_workflow_log_dir`|A directory available to Cromwell|Specifies a path where per-workflow logs will be written. If this is not specified, per-workflow logs will not be copied out of the Cromwell workflow log temporary directory/path before they are deleted.|
|`final_call_logs_dir`|A directory available to Cromwell|Specifies a path where final call logs will be written.  If this is not specified, call logs will not be copied out of the Cromwell workflow execution directory/path.|

Note that these directories should be using the same filesystem as the workflow. Eg if you run on Google's PAPI, you should provide `gs://...` paths.

Example `options.json`:
```json
{
    "final_workflow_outputs_dir": "/Users/michael_scott/cromwell/outputs",
    "use_relative_output_paths": true,
    "final_workflow_log_dir": "/Users/michael_scott/cromwell/wf_logs",
    "final_call_logs_dir": "/Users/michael_scott/cromwell/call_logs"
}
```

With `"use_relative_output_paths": false` (the default) the outputs will look like this

```
final_workflow_outputs_dir/my_workflow/ade68a6d876e8d-8a98d7e9-ad98e9ae8d/call-my_one_task/execution/my_output_picture.jpg
final_workflow_outputs_dir/my_workflow/ade68a6d876e8d-8a98d7e9-ad98e9ae8d/call-my_other_task/execution/created_subdir/submarine.txt
```

The above result will look like this when `"use_relative_output_paths": true`:
```
final_workflow_outputs_dir/my_output_picture.jpg
final_workflow_outputs_dir/created_subdir/submarine.txt
```

This will create file collisions in `final_workflow_outputs_dir` when a workflow is run twice. When cromwell
detects file collisions it will throw an error and report the workflow as failed.

## Call Caching Options

These options can override Cromwell's configured call caching behavior for a single workflow. See the [Call Caching](../CallCaching) section for more details and how to set defaults. The call caching section will also explain how these options interact when, for example, one is set `true` and the other is `false`.

**Note:** If call caching is disabled, these options will be ignored and the options will be treated as though they were `false`.

|Option|Values|Description|
|---|---|---|
|`write_to_cache`|`true` or `false`|If `false`, the completed calls from this workflow will not be added to the cache.  See the [Call Caching](../CallCaching) section for more details.|
|`read_from_cache`|`true` or `false`|If `false`, Cromwell will not search the cache when invoking a call (i.e. every call will be executed unconditionally).  See the [Call Caching](../CallCaching) section for more details.|

Example `options.json`:
```json
{
    "write_to_cache": true,
    "read_from_cache": true
}
```
