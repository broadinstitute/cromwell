# Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi

All URIs are relative to *http://localhost*

| Method | HTTP request | Description |
|------------- | ------------- | -------------|
| [**cancelRun**](Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi.md#cancelRun) | **POST** /api/ga4gh/wes/v1/runs/{run_id}/cancel | Cancel run |
| [**getRunLog**](Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi.md#getRunLog) | **GET** /api/ga4gh/wes/v1/runs/{run_id} | Get run log |
| [**getRunStatus**](Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi.md#getRunStatus) | **GET** /api/ga4gh/wes/v1/runs/{run_id}/status | Get run status |
| [**getServiceInfo**](Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi.md#getServiceInfo) | **GET** /api/ga4gh/wes/v1/service-info | Get service info |
| [**listRuns**](Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi.md#listRuns) | **GET** /api/ga4gh/wes/v1/runs | List runs |
| [**runWorkflow**](Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi.md#runWorkflow) | **POST** /api/ga4gh/wes/v1/runs | Run workflow |


<a name="cancelRun"></a>
# **cancelRun**
> RunId cancelRun(runId)

Cancel run

Cancel a running workflow.

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi apiInstance = new Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi(defaultClient);
    String runId = "runId_example"; // String | Run ID
    try {
      RunId result = apiInstance.cancelRun(runId);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi#cancelRun");
      System.err.println("Status code: " + e.getCode());
      System.err.println("Reason: " + e.getResponseBody());
      System.err.println("Response headers: " + e.getResponseHeaders());
      e.printStackTrace();
    }
  }
}
```

### Parameters

| Name | Type | Description  | Notes |
|------------- | ------------- | ------------- | -------------|
| **runId** | **String**| Run ID | |

### Return type

[**RunId**](RunId.md)

### Authorization

[googleoauth](../README.md#googleoauth)

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

### HTTP response details
| Status code | Description | Response headers |
|-------------|-------------|------------------|
| **200** | Successful Request |  -  |
| **401** | Invalid submission request |  -  |
| **403** | Workflow in terminal status |  -  |
| **404** | Workflow ID Not Found |  -  |
| **500** | Internal Error |  -  |

<a name="getRunLog"></a>
# **getRunLog**
> RunLog getRunLog(runId)

Get run log

This endpoint provides detailed information about a given workflow run. The returned result has information about the outputs produced by this workflow (if available), a log object which allows the stderr and stdout to be retrieved, a log array so stderr/stdout for individual tasks can be retrieved, and the overall state of the workflow run (e.g. RUNNING, see the State section).

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi apiInstance = new Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi(defaultClient);
    String runId = "runId_example"; // String | 
    try {
      RunLog result = apiInstance.getRunLog(runId);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi#getRunLog");
      System.err.println("Status code: " + e.getCode());
      System.err.println("Reason: " + e.getResponseBody());
      System.err.println("Response headers: " + e.getResponseHeaders());
      e.printStackTrace();
    }
  }
}
```

### Parameters

| Name | Type | Description  | Notes |
|------------- | ------------- | ------------- | -------------|
| **runId** | **String**|  | |

### Return type

[**RunLog**](RunLog.md)

### Authorization

[googleoauth](../README.md#googleoauth)

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

### HTTP response details
| Status code | Description | Response headers |
|-------------|-------------|------------------|
| **200** | Successful Request |  -  |
| **401** | Invalid submission request |  -  |
| **403** | Workflow in terminal status |  -  |
| **404** | Workflow ID Not Found |  -  |
| **500** | Internal Error |  -  |

<a name="getRunStatus"></a>
# **getRunStatus**
> RunStatus getRunStatus(runId)

Get run status

This provides an abbreviated (and likely fast depending on implementation) status of the running workflow, returning a simple result with the  overall state of the workflow run (e.g. RUNNING, see the State section).

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi apiInstance = new Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi(defaultClient);
    String runId = "runId_example"; // String | Run ID
    try {
      RunStatus result = apiInstance.getRunStatus(runId);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi#getRunStatus");
      System.err.println("Status code: " + e.getCode());
      System.err.println("Reason: " + e.getResponseBody());
      System.err.println("Response headers: " + e.getResponseHeaders());
      e.printStackTrace();
    }
  }
}
```

### Parameters

| Name | Type | Description  | Notes |
|------------- | ------------- | ------------- | -------------|
| **runId** | **String**| Run ID | |

### Return type

[**RunStatus**](RunStatus.md)

### Authorization

[googleoauth](../README.md#googleoauth)

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

### HTTP response details
| Status code | Description | Response headers |
|-------------|-------------|------------------|
| **200** | Successful Request |  -  |
| **401** | Invalid submission request |  -  |
| **403** | Workflow in terminal status |  -  |
| **404** | Workflow ID Not Found |  -  |
| **500** | Internal Error |  -  |

<a name="getServiceInfo"></a>
# **getServiceInfo**
> ServiceInfo getServiceInfo()

Get service info

May include information related (but not limited to) the workflow descriptor formats, versions supported, the WES API versions supported, and information about general service availability.

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi apiInstance = new Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi(defaultClient);
    try {
      ServiceInfo result = apiInstance.getServiceInfo();
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi#getServiceInfo");
      System.err.println("Status code: " + e.getCode());
      System.err.println("Reason: " + e.getResponseBody());
      System.err.println("Response headers: " + e.getResponseHeaders());
      e.printStackTrace();
    }
  }
}
```

### Parameters
This endpoint does not need any parameter.

### Return type

[**ServiceInfo**](ServiceInfo.md)

### Authorization

[googleoauth](../README.md#googleoauth)

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

### HTTP response details
| Status code | Description | Response headers |
|-------------|-------------|------------------|
| **200** | Successful Request |  -  |
| **500** | Internal Error |  -  |

<a name="listRuns"></a>
# **listRuns**
> RunListResponse listRuns(pageSize, pageToken)

List runs

Runs are listed from newest to oldest. When paging through the list, the client should not make assumptions about live updates, but should assume the contents of the list reflect the workflow list at the moment that the first page is requested. To monitor a specific workflow run, use GetRunStatus or GetRunLog.

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi apiInstance = new Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi(defaultClient);
    Long pageSize = 56L; // Long | OPTIONAL The preferred number of workflow runs to return in a page. If not provided, the implementation should use a default page size. The implementation must not return more items than `page_size`, but it may return fewer.  Clients should not assume that if fewer than `page_size` items are returned that all items have been returned.  The availability of additional pages is indicated by the value of `next_page_token` in the response.
    String pageToken = "pageToken_example"; // String | OPTIONAL Token to use to indicate where to start getting results. If unspecified, return the first page of results.
    try {
      RunListResponse result = apiInstance.listRuns(pageSize, pageToken);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi#listRuns");
      System.err.println("Status code: " + e.getCode());
      System.err.println("Reason: " + e.getResponseBody());
      System.err.println("Response headers: " + e.getResponseHeaders());
      e.printStackTrace();
    }
  }
}
```

### Parameters

| Name | Type | Description  | Notes |
|------------- | ------------- | ------------- | -------------|
| **pageSize** | **Long**| OPTIONAL The preferred number of workflow runs to return in a page. If not provided, the implementation should use a default page size. The implementation must not return more items than &#x60;page_size&#x60;, but it may return fewer.  Clients should not assume that if fewer than &#x60;page_size&#x60; items are returned that all items have been returned.  The availability of additional pages is indicated by the value of &#x60;next_page_token&#x60; in the response. | [optional] |
| **pageToken** | **String**| OPTIONAL Token to use to indicate where to start getting results. If unspecified, return the first page of results. | [optional] |

### Return type

[**RunListResponse**](RunListResponse.md)

### Authorization

[googleoauth](../README.md#googleoauth)

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

### HTTP response details
| Status code | Description | Response headers |
|-------------|-------------|------------------|
| **200** | Successful Request |  -  |
| **400** | Malformed Workflow ID |  -  |
| **401** | Invalid submission request |  -  |
| **403** | Workflow in terminal status |  -  |
| **500** | Internal Error |  -  |

<a name="runWorkflow"></a>
# **runWorkflow**
> RunId runWorkflow(workflowParams, workflowType, workflowTypeVersion, tags, workflowEngineParameters, workflowUrl, workflowAttachment)

Run workflow

This endpoint creates a new workflow run and returns a &#x60;RunId&#x60; to monitor its progress.  The &#x60;workflow_attachment&#x60; array may be used to upload files that are required to execute the workflow, including the primary workflow, tools imported by the workflow, other files referenced by the workflow, or files which are part of the input.  The implementation should stage these files to a temporary directory and execute the workflow from there. These parts must have a Content-Disposition header with a \&quot;filename\&quot; provided for each part.  Filenames may include subdirectories, but must not include references to parent directories with &#39;..&#39; -- implementations should guard against maliciously constructed filenames.  The &#x60;workflow_url&#x60; is either an absolute URL to a workflow file that is accessible by the WES endpoint, or a relative URL corresponding to one of the files attached using &#x60;workflow_attachment&#x60;.  The &#x60;workflow_params&#x60; JSON object specifies input parameters, such as input files.  The exact format of the JSON object depends on the conventions of the workflow language being used.  Input files should either be absolute URLs, or relative URLs corresponding to files uploaded using &#x60;workflow_attachment&#x60;.  The WES endpoint must understand and be able to access URLs supplied in the input.  This is implementation specific.  The &#x60;workflow_type&#x60; is the type of workflow language and must be \&quot;CWL\&quot; or \&quot;WDL\&quot; currently (or another alternative  supported by this WES instance).  The &#x60;workflow_type_version&#x60; is the version of the workflow language submitted and must be one supported by this WES instance.  See the &#x60;RunRequest&#x60; documentation for details about other fields.

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi apiInstance = new Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi(defaultClient);
    String workflowParams = "workflowParams_example"; // String | 
    String workflowType = "workflowType_example"; // String | 
    String workflowTypeVersion = "workflowTypeVersion_example"; // String | 
    String tags = "tags_example"; // String | 
    String workflowEngineParameters = "workflowEngineParameters_example"; // String | 
    String workflowUrl = "workflowUrl_example"; // String | 
    List<File> workflowAttachment = Arrays.asList(); // List<File> | 
    try {
      RunId result = apiInstance.runWorkflow(workflowParams, workflowType, workflowTypeVersion, tags, workflowEngineParameters, workflowUrl, workflowAttachment);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling Ga4GhWorkflowExecutionServiceWesAlphaPreviewApi#runWorkflow");
      System.err.println("Status code: " + e.getCode());
      System.err.println("Reason: " + e.getResponseBody());
      System.err.println("Response headers: " + e.getResponseHeaders());
      e.printStackTrace();
    }
  }
}
```

### Parameters

| Name | Type | Description  | Notes |
|------------- | ------------- | ------------- | -------------|
| **workflowParams** | **String**|  | [optional] |
| **workflowType** | **String**|  | [optional] |
| **workflowTypeVersion** | **String**|  | [optional] |
| **tags** | **String**|  | [optional] |
| **workflowEngineParameters** | **String**|  | [optional] |
| **workflowUrl** | **String**|  | [optional] |
| **workflowAttachment** | **List&lt;File&gt;**|  | [optional] |

### Return type

[**RunId**](RunId.md)

### Authorization

[googleoauth](../README.md#googleoauth)

### HTTP request headers

 - **Content-Type**: multipart/form-data
 - **Accept**: application/json

### HTTP response details
| Status code | Description | Response headers |
|-------------|-------------|------------------|
| **201** | Successful Request |  -  |
| **400** | Malformed Workflow ID |  -  |
| **401** | Invalid submission request |  -  |
| **403** | Workflow in terminal status |  -  |
| **500** | Internal Error |  -  |

