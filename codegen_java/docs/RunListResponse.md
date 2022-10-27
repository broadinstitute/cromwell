

# RunListResponse

The service will return a RunListResponse when receiving a successful RunListRequest.

## Properties

| Name | Type | Description | Notes |
|------------ | ------------- | ------------- | -------------|
|**runs** | [**List&lt;RunStatus&gt;**](RunStatus.md) | A list of workflow runs that the service has executed or is executing. The list is filtered to only include runs that the caller has permission to see. |  [optional] |
|**nextPageToken** | **String** | A token which may be supplied as &#x60;page_token&#x60; in workflow run list request to get the next page of results.  An empty string indicates there are no more items to return. |  [optional] |



