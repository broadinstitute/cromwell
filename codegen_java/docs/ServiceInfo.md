

# ServiceInfo


## Properties

| Name | Type | Description | Notes |
|------------ | ------------- | ------------- | -------------|
|**id** | **String** | Unique ID of this service. Reverse domain name notation is recommended, though not required. The identifier should attempt to be globally unique so it can be used in downstream aggregator services e.g. Service Registry. |  |
|**name** | **String** | Name of this service. Should be human readable. |  |
|**type** | [**ServiceType**](ServiceType.md) |  |  |
|**description** | **String** | Description of the service. Should be human readable and provide information about the service. |  [optional] |
|**organization** | [**ServiceOrganization**](ServiceOrganization.md) |  |  |
|**contactUrl** | **URI** | URL of the contact for the provider of this service, e.g. a link to a contact form (RFC 3986 format), or an email (RFC 2368 format). |  [optional] |
|**documentationUrl** | **URI** | URL of the documentation of this service (RFC 3986 format). This should help someone learn how to use your service, including any specifics required to access data, e.g. authentication. |  [optional] |
|**createdAt** | **OffsetDateTime** | Timestamp describing when the service was first deployed and available (RFC 3339 format) |  [optional] |
|**updatedAt** | **OffsetDateTime** | Timestamp describing when the service was last updated (RFC 3339 format) |  [optional] |
|**environment** | **String** | Environment the service is running in. Use this to distinguish between production, development and testing/staging deployments. Suggested values are prod, test, dev, staging. However this is advised and not enforced. |  [optional] |
|**version** | **String** | Version of the service being described. Semantic versioning is recommended, but other identifiers, such as dates or commit hashes, are also allowed. The version should be changed whenever the service is updated. |  |
|**workflowTypeVersions** | [**Map&lt;String, WorkflowTypeVersion&gt;**](WorkflowTypeVersion.md) |  |  |
|**supportedWesVersions** | **List&lt;String&gt;** | The version(s) of the WES schema supported by this service |  |
|**supportedFilesystemProtocols** | **List&lt;String&gt;** | The filesystem protocols supported by this service, currently these may include common protocols using the terms &#39;http&#39;, &#39;https&#39;, &#39;sftp&#39;, &#39;s3&#39;, &#39;gs&#39;, &#39;file&#39;, or &#39;synapse&#39;, but others  are possible and the terms beyond these core protocols are currently not fixed.   This section reports those protocols (either common or not) supported by this WES service. |  |
|**workflowEngineVersions** | **Map&lt;String, String&gt;** |  |  |
|**defaultWorkflowEngineParameters** | [**List&lt;DefaultWorkflowEngineParameter&gt;**](DefaultWorkflowEngineParameter.md) | Each workflow engine can present additional parameters that can be sent to the workflow engine. This message will list the default values, and their types for each workflow engine. |  |
|**systemStateCounts** | **Map&lt;String, Long&gt;** |  |  |
|**authInstructionsUrl** | **String** | A web page URL with human-readable instructions on how to get an authorization token for use with a specific WES endpoint. |  |
|**tags** | **Map&lt;String, String&gt;** |  |  |



