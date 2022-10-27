

# Log

Log and other info

## Properties

| Name | Type | Description | Notes |
|------------ | ------------- | ------------- | -------------|
|**name** | **String** | The task or workflow name |  [optional] |
|**cmd** | **List&lt;String&gt;** | The command line that was executed |  [optional] |
|**startTime** | **String** | When the command started executing, in ISO 8601 format \&quot;%Y-%m-%dT%H:%M:%SZ\&quot; |  [optional] |
|**endTime** | **String** | When the command stopped executing (completed, failed, or cancelled), in ISO 8601 format \&quot;%Y-%m-%dT%H:%M:%SZ\&quot; |  [optional] |
|**stdout** | **String** | A URL to retrieve standard output logs of the workflow run or task.  This URL may change between status requests, or may not be available until the task or workflow has finished execution.  Should be available using the same credentials used to access the WES endpoint. |  [optional] |
|**stderr** | **String** | A URL to retrieve standard error logs of the workflow run or task.  This URL may change between status requests, or may not be available until the task or workflow has finished execution.  Should be available using the same credentials used to access the WES endpoint. |  [optional] |
|**exitCode** | **Integer** | Exit code of the program |  [optional] |



