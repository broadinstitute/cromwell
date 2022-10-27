

# CallMetadata

Call level metadata

## Properties

| Name | Type | Description | Notes |
|------------ | ------------- | ------------- | -------------|
|**inputs** | **Object** | Mapping of input fully qualified names to stringified values |  |
|**executionStatus** | **String** | Status in Cromwell execution terms. |  |
|**backend** | **String** | The type of backend on which the call executed (e.g. JES, SGE, Local) |  [optional] |
|**backendStatus** | **String** | Status in backend-specific terms.  Currently this will only be defined for the JES backend. |  [optional] |
|**start** | **OffsetDateTime** | Start datetime of the call execution in ISO8601 format with milliseconds |  [optional] |
|**end** | **OffsetDateTime** | End datetime of the call execution in ISO8601 format with milliseconds |  [optional] |
|**jobId** | **String** | Backend-specific job ID |  [optional] |
|**failures** | [**FailureMessage**](FailureMessage.md) |  |  [optional] |
|**returnCode** | **Integer** | Call execution return code |  [optional] |
|**stdout** | **String** | Path to the standard output file for this call |  [optional] |
|**stderr** | **String** | Path to the standard error file for this call |  [optional] |
|**backendLogs** | **Object** | Paths to backend specific logs for this call |  [optional] |



