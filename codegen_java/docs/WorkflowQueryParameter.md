

# WorkflowQueryParameter

Workflow query parameters

## Properties

| Name | Type | Description | Notes |
|------------ | ------------- | ------------- | -------------|
|**submission** | **OffsetDateTime** | Returns only workflows with an equal or later submission time. Can be specified at most once. If both submission time and start date are specified, submission time should be before or equal to start date.  |  [optional] |
|**start** | **OffsetDateTime** | Returns only workflows with an equal or later start datetime.  Can be specified at most once. If both start and end date are specified, start date must be before or equal to end date.  |  [optional] |
|**end** | **OffsetDateTime** | Returns only workflows with an equal or earlier end datetime.  Can be specified at most once. If both start and end date are specified, start date must be before or equal to end date.  |  [optional] |
|**status** | [**StatusEnum**](#StatusEnum) | Returns only workflows with the specified status.  If specified multiple times, returns workflows in any of the specified statuses.  |  [optional] |
|**name** | **String** | Returns only workflows with the specified name.  If specified multiple times, returns workflows with any of the specified names.  |  [optional] |
|**id** | **String** | Returns only workflows with the specified workflow id.  If specified multiple times, returns workflows with any of the specified workflow ids.  |  [optional] |
|**excludeLabelAnd** | **List** | Excludes workflows with the specified label.  If specified multiple times, excludes workflows with all of the specified label keys. Specify the label key and label value pair as separated with \&quot;label-key:label-value\&quot;  |  [optional] |
|**excludeLabelOr** | **List** | Excludes workflows with the specified label.  If specified multiple times, excludes workflows with any of the specified label keys. Specify the label key and label value pair as separated with \&quot;label-key:label-value\&quot;  |  [optional] |
|**includeSubworkflows** | **Boolean** | Include subworkflows in results. By default, it is taken as true. |  [optional] |
|**page** | **Integer** | When pageSize is set, which page of results to return. If not set, the first page of &#39;pageSize&#39; results will be returned. |  [optional] |
|**pageSize** | **Integer** | The number of results to return at a time |  [optional] |



## Enum: StatusEnum

| Name | Value |
|---- | -----|
| SUBMITTED | &quot;Submitted&quot; |
| RUNNING | &quot;Running&quot; |
| ABORTING | &quot;Aborting&quot; |
| FAILED | &quot;Failed&quot; |
| SUCCEEDED | &quot;Succeeded&quot; |
| ABORTED | &quot;Aborted&quot; |



