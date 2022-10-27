

# WorkflowMetadataResponse

Workflow and call level metadata

## Properties

| Name | Type | Description | Notes |
|------------ | ------------- | ------------- | -------------|
|**id** | **String** | The identifier of the workflow |  |
|**status** | **String** | The status of the workflow |  |
|**submission** | **OffsetDateTime** | Submission datetime of the workflow in ISO8601 format with milliseconds |  |
|**start** | **OffsetDateTime** | Start datetime of the workflow in ISO8601 format with milliseconds |  [optional] |
|**end** | **OffsetDateTime** | End datetime of the workflow in ISO8601 format with milliseconds |  [optional] |
|**inputs** | **Object** | Map of input keys to input values |  [optional] |
|**outputs** | **Object** | Map of output keys to output values |  [optional] |
|**calls** | [**CallMetadata**](CallMetadata.md) |  |  [optional] |
|**failures** | [**FailureMessage**](FailureMessage.md) |  |  [optional] |



