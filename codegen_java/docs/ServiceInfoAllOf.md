

# ServiceInfoAllOf


## Properties

| Name | Type | Description | Notes |
|------------ | ------------- | ------------- | -------------|
|**workflowTypeVersions** | [**Map&lt;String, WorkflowTypeVersion&gt;**](WorkflowTypeVersion.md) |  |  |
|**supportedWesVersions** | **List&lt;String&gt;** | The version(s) of the WES schema supported by this service |  |
|**supportedFilesystemProtocols** | **List&lt;String&gt;** | The filesystem protocols supported by this service, currently these may include common protocols using the terms &#39;http&#39;, &#39;https&#39;, &#39;sftp&#39;, &#39;s3&#39;, &#39;gs&#39;, &#39;file&#39;, or &#39;synapse&#39;, but others  are possible and the terms beyond these core protocols are currently not fixed.   This section reports those protocols (either common or not) supported by this WES service. |  |
|**workflowEngineVersions** | **Map&lt;String, String&gt;** |  |  |
|**defaultWorkflowEngineParameters** | [**List&lt;DefaultWorkflowEngineParameter&gt;**](DefaultWorkflowEngineParameter.md) | Each workflow engine can present additional parameters that can be sent to the workflow engine. This message will list the default values, and their types for each workflow engine. |  |
|**systemStateCounts** | **Map&lt;String, Long&gt;** |  |  |
|**authInstructionsUrl** | **String** | A web page URL with human-readable instructions on how to get an authorization token for use with a specific WES endpoint. |  |
|**tags** | **Map&lt;String, String&gt;** |  |  |



