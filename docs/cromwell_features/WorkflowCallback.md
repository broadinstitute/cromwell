The workflow callback is a simple way to integrate Cromwell with an external system. When each workflow reaches a terminal
state, Cromwell will attempt to POST a message to a provided URL (see below for schema of this message). 
Messages are sent for root workflows only, not subworkflows. Callback status information, including success or failure, 
will be recorded in workflow metadata with keys containing `workflowCallback`.

### Configuration

This feature will only be used if enabled via config. All config items except `enabled` are optional.

```
workflow-state-callback {
  enabled: true
  num-threads: 5
  endpoint: "http://example.com"
  auth.azure: true
  request-backoff {
    min: "3 seconds",
    max: "5 minutes",
    multiplier: 1.1
  }
  max-retries = 10
}
```

 * `enabled`: This boolean controls whether a callback will be attempted or not.
 * `num-threads`: The number of threads Cromwell will allocate for performing callbacks.
 * `endpoint`: This is the default URL to send the message to. If this is unset, and no URL is set in workflow options, no callback will be sent.
 * `auth.azure`: If true, and if Cromwell is running in an Azure environment, Cromwell will include an auth header with bearer token generated from local Azure credentials.
 * `request-backoff` and `max-retries`: Include these to override the default retry behavior (default behavior shown here).

### Workflow Options

You may choose to override the `endpoint` set in config by including this workflow option:
```
{
  "workflow_callback_uri": "http://mywebsite.com"
}
```

### Callback schema

Below is an example of a callback request body.

```
{
  "workflowId": "00001111-2222-3333-4444-555566667777",
  "state": "Succeeded",
  "outputs": {
    "task1.out": 5,
    "task2.out": "/some/file.txt"
  }
}
```

 * `workflowId`: The UUID of the workflow
 * `state`: The terminal state of the workflow. The list of possible values is: `Succeeded`, `Failed`, `Aborted`
 * `outputs`: The final outputs of the workflow, as would be returned from the `api/workflows/{version}/{id}/outputs` endpoint. Expected to be empty when the workflow is not successful..
 * `failures`: A list of strings describing the workflow's failures. Expected to be empty if the workflow did not fail.
