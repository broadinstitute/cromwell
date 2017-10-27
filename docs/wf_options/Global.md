# Global Workflow Options 

The following workflow options apply to all workflows and their calls regardless of the backend being used.
For information on how to compose these into a workflow options file see [Workflow Options](Overview)

## Runtime Attributes

Some options allow you to override or set defaults for runtime attributes.

### Setting Default Runtime Attributes

You can supply a default for any [Runtime Attributes](RuntimeAttributes) by adding a `default_runtime_attributes` map to your workflow options file. Use the key to provide the attribute name and the value to supply the default. 

These defaults replace any defaults in the Cromwell configuration file but are themselves replaced by any values explicitly provided by the task in the WDL or CWL file.

Example:
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

Example:
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

Example:
```json
{
    "workflow_failure_mode": "ContinueWhilePossible"
}
```

## Output Copying
|Option|Value|Description|
|---|---|---|
|`final_workflow_outputs_dir`|A directory available to Cromwell|Specifies a path where final workflow outputs will be written. If this is not specified, workflow outputs will not be copied out of the Cromwell workflow execution directory/path.|
|`final_workflow_log_dir`|A directory available to Cromwell|Specifies a path where per-workflow logs will be written. If this is not specified, per-workflow logs will not be copied out of the Cromwell workflow log temporary directory/path before they are deleted. **TODO**: are they actually deleted?!?|
|`final_call_logs_dir`|A directory available to Cromwell|Specifies a path where final call logs will be written.  If this is not specified, call logs will not be copied out of the Cromwell workflow execution directory/path.|

Note that these directories should be using the same filesystem as the workflow. Eg if you run on Google's PAPI, you should provide `gs://...` paths.

Example:
```json
{
    "final_workflow_outputs_dir": "/Users/michael_scott/cromwell/outputs",
    "final_workflow_log_dir": "/Users/michael_scott/cromwell/wf_logs",
    "final_call_logs_dir": "/Users/michael_scott/cromwell/call_logs"
}
```

## Call Caching Options

These options can override Cromwell's configured call caching behavior for a single workflow. See the [Call Caching](CallCaching) section for more details and how to set defaults. The call caching section will also explain how these options interact when, for example, one is set `true` and the other is `false`.

|Option|Values|Description|
|---|---|---|
|`write_to_cache`|`true` or `false`|If `false`, the completed calls from this workflow will not be added to the cache.  See the [Call Caching](CallCaching) section for more details.|
|`read_from_cache`|`true` or `false`|If `false`, Cromwell will not search the cache when invoking a call (i.e. every call will be executed unconditionally).  See the [Call Caching](CallCaching) section for more details.|

Example:
```json
{
    "write_to_cache": true,
    "read_from_cache": true
}
```

