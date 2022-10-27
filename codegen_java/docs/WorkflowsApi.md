# WorkflowsApi

All URIs are relative to *http://localhost*

| Method | HTTP request | Description |
|------------- | ------------- | -------------|
| [**abort**](WorkflowsApi.md#abort) | **POST** /api/workflows/{version}/{id}/abort | Abort a running workflow |
| [**backends**](WorkflowsApi.md#backends) | **GET** /api/workflows/{version}/backends | List the supported backends |
| [**callCacheDiff**](WorkflowsApi.md#callCacheDiff) | **GET** /api/workflows/{version}/callcaching/diff | Explain hashing differences for 2 calls |
| [**labels**](WorkflowsApi.md#labels) | **GET** /api/workflows/{version}/{id}/labels | Retrieves the current labels for a workflow |
| [**logs**](WorkflowsApi.md#logs) | **GET** /api/workflows/{version}/{id}/logs | Get the logs for a workflow |
| [**metadata**](WorkflowsApi.md#metadata) | **GET** /api/workflows/{version}/{id}/metadata | Get workflow and call-level metadata for a specified workflow |
| [**outputs**](WorkflowsApi.md#outputs) | **GET** /api/workflows/{version}/{id}/outputs | Get the outputs for a workflow |
| [**queryGet**](WorkflowsApi.md#queryGet) | **GET** /api/workflows/{version}/query | Get workflows matching some criteria |
| [**queryPost**](WorkflowsApi.md#queryPost) | **POST** /api/workflows/{version}/query | Get workflows matching some criteria |
| [**releaseHold**](WorkflowsApi.md#releaseHold) | **POST** /api/workflows/{version}/{id}/releaseHold | Switch a workflow from &#39;On Hold&#39; to &#39;Submitted&#39; status |
| [**status**](WorkflowsApi.md#status) | **GET** /api/workflows/{version}/{id}/status | Retrieves the current state for a workflow |
| [**submit**](WorkflowsApi.md#submit) | **POST** /api/workflows/{version} | Submit a workflow for execution |
| [**submitBatch**](WorkflowsApi.md#submitBatch) | **POST** /api/workflows/{version}/batch | Submit a batch of workflows for execution |
| [**timing**](WorkflowsApi.md#timing) | **GET** /api/workflows/{version}/{id}/timing | Get a visual diagram of a running workflow |
| [**updateLabels**](WorkflowsApi.md#updateLabels) | **PATCH** /api/workflows/{version}/{id}/labels | Update labels for a workflow |


<a name="abort"></a>
# **abort**
> WorkflowIdAndStatus abort(version, id)

Abort a running workflow

Request Cromwell to abort a running workflow. For instance this might be necessary in cases where you have submitted a workflow with incorrect inputs or no longer need the results. Cromwell will schedule a halt of all currently running jobs from this workflow.

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.WorkflowsApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    WorkflowsApi apiInstance = new WorkflowsApi(defaultClient);
    String version = "v1"; // String | Cromwell API Version
    String id = "id_example"; // String | A workflow ID
    try {
      WorkflowIdAndStatus result = apiInstance.abort(version, id);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling WorkflowsApi#abort");
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
| **version** | **String**| Cromwell API Version | [default to v1] |
| **id** | **String**| A workflow ID | |

### Return type

[**WorkflowIdAndStatus**](WorkflowIdAndStatus.md)

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
| **403** | Workflow in terminal status |  -  |
| **404** | Workflow ID Not Found |  -  |
| **500** | Internal Error |  -  |

<a name="backends"></a>
# **backends**
> BackendResponse backends(version)

List the supported backends

Returns the backends supported by this Cromwell server, as well as the default backend.

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.WorkflowsApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    WorkflowsApi apiInstance = new WorkflowsApi(defaultClient);
    String version = "v1"; // String | Cromwell API Version
    try {
      BackendResponse result = apiInstance.backends(version);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling WorkflowsApi#backends");
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
| **version** | **String**| Cromwell API Version | [default to v1] |

### Return type

[**BackendResponse**](BackendResponse.md)

### Authorization

[googleoauth](../README.md#googleoauth)

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

### HTTP response details
| Status code | Description | Response headers |
|-------------|-------------|------------------|
| **200** | Successful Request |  -  |

<a name="callCacheDiff"></a>
# **callCacheDiff**
> WorkflowIdAndStatus callCacheDiff(version, workflowA, callA, workflowB, callB, indexA, indexB)

Explain hashing differences for 2 calls

This endpoint returns the hash differences between 2 completed (successfully or not) calls.

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.WorkflowsApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    WorkflowsApi apiInstance = new WorkflowsApi(defaultClient);
    String version = "v1"; // String | Cromwell API Version
    String workflowA = "workflowA_example"; // String | Workflow Id of the first workflow
    String callA = "callA_example"; // String | Fully qualified name, including workflow name, of the first call.
    String workflowB = "workflowB_example"; // String | Workflow Id of the second workflow
    String callB = "callB_example"; // String | Fully qualified name, including workflow name, of the second call
    Integer indexA = 56; // Integer | Shard index for the first call for cases where the requested call was part of a scatter.
    Integer indexB = 56; // Integer | Shard index for the second call for cases where the requested call was part of a scatter.
    try {
      WorkflowIdAndStatus result = apiInstance.callCacheDiff(version, workflowA, callA, workflowB, callB, indexA, indexB);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling WorkflowsApi#callCacheDiff");
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
| **version** | **String**| Cromwell API Version | [default to v1] |
| **workflowA** | **String**| Workflow Id of the first workflow | |
| **callA** | **String**| Fully qualified name, including workflow name, of the first call. | |
| **workflowB** | **String**| Workflow Id of the second workflow | |
| **callB** | **String**| Fully qualified name, including workflow name, of the second call | |
| **indexA** | **Integer**| Shard index for the first call for cases where the requested call was part of a scatter. | [optional] |
| **indexB** | **Integer**| Shard index for the second call for cases where the requested call was part of a scatter. | [optional] |

### Return type

[**WorkflowIdAndStatus**](WorkflowIdAndStatus.md)

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
| **404** | No matching cache entry. Cromwell versions prior to 28 will not have recorded information necessary for this endpoint and thus will also appear to not exist. |  -  |
| **500** | Internal Error |  -  |

<a name="labels"></a>
# **labels**
> LabelsResponse labels(version, id)

Retrieves the current labels for a workflow

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.WorkflowsApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    WorkflowsApi apiInstance = new WorkflowsApi(defaultClient);
    String version = "v1"; // String | Cromwell API Version
    String id = "id_example"; // String | A workflow ID
    try {
      LabelsResponse result = apiInstance.labels(version, id);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling WorkflowsApi#labels");
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
| **version** | **String**| Cromwell API Version | [default to v1] |
| **id** | **String**| A workflow ID | |

### Return type

[**LabelsResponse**](LabelsResponse.md)

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
| **404** | Workflow ID Not Found |  -  |
| **500** | Internal Error |  -  |

<a name="logs"></a>
# **logs**
> WorkflowIdAndStatus logs(version, id)

Get the logs for a workflow

Returns paths to the standard out and standard error files that were generated during the execution of all calls in a workflow. A call has one or more standard out and standard error logs, depending on if the call was scattered or not. In the latter case, one log is provided for each instance of the call that has been run.

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.WorkflowsApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    WorkflowsApi apiInstance = new WorkflowsApi(defaultClient);
    String version = "v1"; // String | Cromwell API Version
    String id = "id_example"; // String | A workflow ID
    try {
      WorkflowIdAndStatus result = apiInstance.logs(version, id);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling WorkflowsApi#logs");
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
| **version** | **String**| Cromwell API Version | [default to v1] |
| **id** | **String**| A workflow ID | |

### Return type

[**WorkflowIdAndStatus**](WorkflowIdAndStatus.md)

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
| **404** | Workflow ID Not Found |  -  |
| **500** | Internal Error |  -  |

<a name="metadata"></a>
# **metadata**
> WorkflowMetadataResponse metadata(version, id, includeKey, excludeKey, expandSubWorkflows)

Get workflow and call-level metadata for a specified workflow

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.WorkflowsApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    WorkflowsApi apiInstance = new WorkflowsApi(defaultClient);
    String version = "v1"; // String | Cromwell API Version
    String id = "id_example"; // String | A workflow ID
    List<String> includeKey = Arrays.asList(); // List<String> | When specified, filters metadata to only return fields with names which begins with this value. This key is used relative to the root of the response *and* relative to each call's metadata fields. 
    List<String> excludeKey = Arrays.asList(); // List<String> | When specified, filters metadata to not return any field with a name which begins with this value. This key is used relative to the root of the response *and* relative to each call's metadata fields. Use 'calls' to filter out all call level metadata. 
    Boolean expandSubWorkflows = true; // Boolean | When true, metadata for sub workflows will be fetched and inserted automatically in the metadata response. 
    try {
      WorkflowMetadataResponse result = apiInstance.metadata(version, id, includeKey, excludeKey, expandSubWorkflows);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling WorkflowsApi#metadata");
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
| **version** | **String**| Cromwell API Version | [default to v1] |
| **id** | **String**| A workflow ID | |
| **includeKey** | [**List&lt;String&gt;**](String.md)| When specified, filters metadata to only return fields with names which begins with this value. This key is used relative to the root of the response *and* relative to each call&#39;s metadata fields.  | [optional] |
| **excludeKey** | [**List&lt;String&gt;**](String.md)| When specified, filters metadata to not return any field with a name which begins with this value. This key is used relative to the root of the response *and* relative to each call&#39;s metadata fields. Use &#39;calls&#39; to filter out all call level metadata.  | [optional] |
| **expandSubWorkflows** | **Boolean**| When true, metadata for sub workflows will be fetched and inserted automatically in the metadata response.  | [optional] |

### Return type

[**WorkflowMetadataResponse**](WorkflowMetadataResponse.md)

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
| **404** | Workflow ID Not Found |  -  |
| **500** | Internal Error |  -  |

<a name="outputs"></a>
# **outputs**
> WorkflowIdAndStatus outputs(version, id)

Get the outputs for a workflow

Retrieve the outputs for the specified workflow. Cromwell will return any outputs which currently exist even if a workflow has not successfully completed.

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.WorkflowsApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    WorkflowsApi apiInstance = new WorkflowsApi(defaultClient);
    String version = "v1"; // String | Cromwell API Version
    String id = "id_example"; // String | A workflow ID
    try {
      WorkflowIdAndStatus result = apiInstance.outputs(version, id);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling WorkflowsApi#outputs");
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
| **version** | **String**| Cromwell API Version | [default to v1] |
| **id** | **String**| A workflow ID | |

### Return type

[**WorkflowIdAndStatus**](WorkflowIdAndStatus.md)

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
| **404** | Workflow ID Not Found |  -  |
| **500** | Internal Error |  -  |

<a name="queryGet"></a>
# **queryGet**
> WorkflowQueryResponse queryGet(version, submission, start, end, status, name, id, label, labelor, excludeLabelAnd, excludeLabelOr, additionalQueryResultFields, includeSubworkflows)

Get workflows matching some criteria

Query for workflows which match various criteria. When a combination of criteria are applied the endpoint will return

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.WorkflowsApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    WorkflowsApi apiInstance = new WorkflowsApi(defaultClient);
    String version = "v1"; // String | Cromwell API Version
    OffsetDateTime submission = OffsetDateTime.now(); // OffsetDateTime | Returns only workflows with an equal or later submission time. Can be specified at most once. If both submission time and start date are specified, submission time should be before or equal to start date. 
    OffsetDateTime start = OffsetDateTime.now(); // OffsetDateTime | Returns only workflows with an equal or later start datetime.  Can be specified at most once. If both start and end date are specified, start date must be before or equal to end date. 
    OffsetDateTime end = OffsetDateTime.now(); // OffsetDateTime | Returns only workflows with an equal or earlier end datetime.  Can be specified at most once. If both start and end date are specified, start date must be before or equal to end date. 
    List<String> status = Arrays.asList(); // List<String> | Returns only workflows with the specified status.  If specified multiple times, returns workflows in any of the specified statuses. 
    List<String> name = Arrays.asList(); // List<String> | Returns only workflows with the specified name.  If specified multiple times, returns workflows with any of the specified names. 
    List<String> id = Arrays.asList(); // List<String> | Returns only workflows with the specified workflow id.  If specified multiple times, returns workflows with any of the specified workflow ids. 
    List<String> label = Arrays.asList(); // List<String> | Returns workflows with the specified label keys.  If specified multiple times, returns workflows with all of the specified label keys. Specify the label key and label value pair as separated with \"label-key:label-value\" 
    List<String> labelor = Arrays.asList(); // List<String> | Returns workflows with the specified label keys.  If specified multiple times, returns workflows with any of the specified label keys. Specify the label key and label value pair as separated with \"label-key:label-value\" 
    List<String> excludeLabelAnd = Arrays.asList(); // List<String> | Excludes workflows with the specified label.  If specified multiple times, excludes workflows with all of the specified label keys. Specify the label key and label value pair as separated with \"label-key:label-value\" 
    List<String> excludeLabelOr = Arrays.asList(); // List<String> | Excludes workflows with the specified label.  If specified multiple times, excludes workflows with any of the specified label keys. Specify the label key and label value pair as separated with \"label-key:label-value\" 
    List<String> additionalQueryResultFields = Arrays.asList(); // List<String> | Currently only 'labels' is a valid value here. Use it to include a list of labels with each result. 
    Boolean includeSubworkflows = true; // Boolean | Include subworkflows in results. By default, it is taken as true.
    try {
      WorkflowQueryResponse result = apiInstance.queryGet(version, submission, start, end, status, name, id, label, labelor, excludeLabelAnd, excludeLabelOr, additionalQueryResultFields, includeSubworkflows);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling WorkflowsApi#queryGet");
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
| **version** | **String**| Cromwell API Version | [default to v1] |
| **submission** | **OffsetDateTime**| Returns only workflows with an equal or later submission time. Can be specified at most once. If both submission time and start date are specified, submission time should be before or equal to start date.  | [optional] |
| **start** | **OffsetDateTime**| Returns only workflows with an equal or later start datetime.  Can be specified at most once. If both start and end date are specified, start date must be before or equal to end date.  | [optional] |
| **end** | **OffsetDateTime**| Returns only workflows with an equal or earlier end datetime.  Can be specified at most once. If both start and end date are specified, start date must be before or equal to end date.  | [optional] |
| **status** | [**List&lt;String&gt;**](String.md)| Returns only workflows with the specified status.  If specified multiple times, returns workflows in any of the specified statuses.  | [optional] [enum: Submitted, Running, Aborting, Failed, Succeeded, Aborted] |
| **name** | [**List&lt;String&gt;**](String.md)| Returns only workflows with the specified name.  If specified multiple times, returns workflows with any of the specified names.  | [optional] |
| **id** | [**List&lt;String&gt;**](String.md)| Returns only workflows with the specified workflow id.  If specified multiple times, returns workflows with any of the specified workflow ids.  | [optional] |
| **label** | [**List&lt;String&gt;**](String.md)| Returns workflows with the specified label keys.  If specified multiple times, returns workflows with all of the specified label keys. Specify the label key and label value pair as separated with \&quot;label-key:label-value\&quot;  | [optional] |
| **labelor** | [**List&lt;String&gt;**](String.md)| Returns workflows with the specified label keys.  If specified multiple times, returns workflows with any of the specified label keys. Specify the label key and label value pair as separated with \&quot;label-key:label-value\&quot;  | [optional] |
| **excludeLabelAnd** | [**List&lt;String&gt;**](String.md)| Excludes workflows with the specified label.  If specified multiple times, excludes workflows with all of the specified label keys. Specify the label key and label value pair as separated with \&quot;label-key:label-value\&quot;  | [optional] |
| **excludeLabelOr** | [**List&lt;String&gt;**](String.md)| Excludes workflows with the specified label.  If specified multiple times, excludes workflows with any of the specified label keys. Specify the label key and label value pair as separated with \&quot;label-key:label-value\&quot;  | [optional] |
| **additionalQueryResultFields** | [**List&lt;String&gt;**](String.md)| Currently only &#39;labels&#39; is a valid value here. Use it to include a list of labels with each result.  | [optional] |
| **includeSubworkflows** | **Boolean**| Include subworkflows in results. By default, it is taken as true. | [optional] |

### Return type

[**WorkflowQueryResponse**](WorkflowQueryResponse.md)

### Authorization

[googleoauth](../README.md#googleoauth)

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

### HTTP response details
| Status code | Description | Response headers |
|-------------|-------------|------------------|
| **200** | Successful Request |  -  |
| **403** | Workflow in terminal status |  -  |
| **404** | Workflow ID Not Found |  -  |
| **500** | Internal Error |  -  |

<a name="queryPost"></a>
# **queryPost**
> WorkflowQueryResponse queryPost(version, parameters)

Get workflows matching some criteria

Query workflows by start dates, end dates, names, ids, labels, or statuses.

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.WorkflowsApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    WorkflowsApi apiInstance = new WorkflowsApi(defaultClient);
    String version = "v1"; // String | Cromwell API Version
    List<WorkflowQueryParameter> parameters = Arrays.asList(); // List<WorkflowQueryParameter> | Same query parameters as GET /query endpoint, submitted as a json list. Example: [{\"status\":\"Success\"},{\"status\":\"Failed\"}] 
    try {
      WorkflowQueryResponse result = apiInstance.queryPost(version, parameters);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling WorkflowsApi#queryPost");
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
| **version** | **String**| Cromwell API Version | [default to v1] |
| **parameters** | [**List&lt;WorkflowQueryParameter&gt;**](WorkflowQueryParameter.md)| Same query parameters as GET /query endpoint, submitted as a json list. Example: [{\&quot;status\&quot;:\&quot;Success\&quot;},{\&quot;status\&quot;:\&quot;Failed\&quot;}]  | |

### Return type

[**WorkflowQueryResponse**](WorkflowQueryResponse.md)

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
| **500** | Internal Error |  -  |

<a name="releaseHold"></a>
# **releaseHold**
> WorkflowIdAndStatus releaseHold(version, id)

Switch a workflow from &#39;On Hold&#39; to &#39;Submitted&#39; status

Request Cromwell to release the hold on a workflow. It will switch the status of a workflow from &#39;On Hold&#39; to &#39;Submitted&#39; so it can be picked for running. For instance this might be necessary in cases where you have submitted a workflow with workflowOnHold &#x3D; true.

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.WorkflowsApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    WorkflowsApi apiInstance = new WorkflowsApi(defaultClient);
    String version = "v1"; // String | Cromwell API Version
    String id = "id_example"; // String | A workflow ID
    try {
      WorkflowIdAndStatus result = apiInstance.releaseHold(version, id);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling WorkflowsApi#releaseHold");
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
| **version** | **String**| Cromwell API Version | [default to v1] |
| **id** | **String**| A workflow ID | |

### Return type

[**WorkflowIdAndStatus**](WorkflowIdAndStatus.md)

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
| **403** | Workflow in terminal status |  -  |
| **404** | Workflow ID Not Found |  -  |
| **500** | Internal Error |  -  |

<a name="status"></a>
# **status**
> WorkflowIdAndStatus status(version, id)

Retrieves the current state for a workflow

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.WorkflowsApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    WorkflowsApi apiInstance = new WorkflowsApi(defaultClient);
    String version = "v1"; // String | Cromwell API Version
    String id = "id_example"; // String | A workflow ID
    try {
      WorkflowIdAndStatus result = apiInstance.status(version, id);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling WorkflowsApi#status");
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
| **version** | **String**| Cromwell API Version | [default to v1] |
| **id** | **String**| A workflow ID | |

### Return type

[**WorkflowIdAndStatus**](WorkflowIdAndStatus.md)

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
| **404** | Workflow ID Not Found |  -  |
| **500** | Internal Error |  -  |

<a name="submit"></a>
# **submit**
> WorkflowIdAndStatus submit(version, workflowSource, workflowUrl, workflowOnHold, workflowInputs, workflowInputs2, workflowInputs3, workflowInputs4, workflowInputs5, workflowOptions, workflowType, workflowRoot, workflowTypeVersion, labels, workflowDependencies, requestedWorkflowId)

Submit a workflow for execution

Submits a workflow to Cromwell. Note that this endpoint can accept an unlimited number of input files via workflowInputs_N but swagger needs them to be explicitly defined so we have provided 5 as an example.

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.WorkflowsApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    WorkflowsApi apiInstance = new WorkflowsApi(defaultClient);
    String version = "v1"; // String | Cromwell API Version
    String workflowSource = "workflowSource_example"; // String | The workflow source file to submit for execution. Either workflow source or workflow url is required.
    String workflowUrl = "workflowUrl_example"; // String | URL which points to the workflow. Either workflow source or workflow url is required.
    Boolean workflowOnHold = true; // Boolean | Put workflow on hold upon submission. By default, it is taken as false.
    String workflowInputs = "workflowInputs_example"; // String | JSON or YAML file containing the inputs as an object. For WDL workflows a skeleton file can be generated from WOMtool using the \\\"inputs\\\" subcommand. When multiple files are specified, in case of key conflicts between multiple input JSON files, higher values of x in workflowInputs_x override lower values. For example, an input specified in workflowInputs_3 will override an input with the same name in workflowInputs or workflowInputs_2. Similarly, an input key specified in workflowInputs_5 will override an identical input key in any other input file.
    String workflowInputs2 = "workflowInputs2_example"; // String | A second JSON or YAML file containing inputs.
    String workflowInputs3 = "workflowInputs3_example"; // String | A third JSON or YAML file containing inputs.
    String workflowInputs4 = "workflowInputs4_example"; // String | A fourth JSON or YAML file containing inputs.
    String workflowInputs5 = "workflowInputs5_example"; // String | A fifth JSON or YAML file containing inputs.
    String workflowOptions = "workflowOptions_example"; // String | JSON file containing configuration options for the execution of this workflow.
    String workflowType = "WDL"; // String | The workflow language for the file you submitted. Cromwell currently supports WDL and CWL.
    String workflowRoot = "workflowRoot_example"; // String | The root object to be run. Only necessary for CWL submissions containing multiple objects (in an array).
    String workflowTypeVersion = "draft-2"; // String | The specification version for the workflow language being used. For WDL, Cromwell currently supports draft-2 and 1.0. For CWL, Cromwell currently supports v1.0.
    String labels = "labels_example"; // String | JSON object of labels to apply to this workflow.
    String workflowDependencies = "workflowDependencies_example"; // String | ZIP file containing workflow source files that are used to resolve local imports. This zip bundle will be unpacked in a sandbox accessible to this workflow.
    String requestedWorkflowId = "requestedWorkflowId_example"; // String | An ID to ascribe to this workflow. Must be a JSON string in UUID-format. If not supplied a random ID will be generated for the workflow.
    try {
      WorkflowIdAndStatus result = apiInstance.submit(version, workflowSource, workflowUrl, workflowOnHold, workflowInputs, workflowInputs2, workflowInputs3, workflowInputs4, workflowInputs5, workflowOptions, workflowType, workflowRoot, workflowTypeVersion, labels, workflowDependencies, requestedWorkflowId);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling WorkflowsApi#submit");
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
| **version** | **String**| Cromwell API Version | [default to v1] |
| **workflowSource** | **String**| The workflow source file to submit for execution. Either workflow source or workflow url is required. | [optional] |
| **workflowUrl** | **String**| URL which points to the workflow. Either workflow source or workflow url is required. | [optional] |
| **workflowOnHold** | **Boolean**| Put workflow on hold upon submission. By default, it is taken as false. | [optional] |
| **workflowInputs** | **String**| JSON or YAML file containing the inputs as an object. For WDL workflows a skeleton file can be generated from WOMtool using the \\\&quot;inputs\\\&quot; subcommand. When multiple files are specified, in case of key conflicts between multiple input JSON files, higher values of x in workflowInputs_x override lower values. For example, an input specified in workflowInputs_3 will override an input with the same name in workflowInputs or workflowInputs_2. Similarly, an input key specified in workflowInputs_5 will override an identical input key in any other input file. | [optional] |
| **workflowInputs2** | **String**| A second JSON or YAML file containing inputs. | [optional] |
| **workflowInputs3** | **String**| A third JSON or YAML file containing inputs. | [optional] |
| **workflowInputs4** | **String**| A fourth JSON or YAML file containing inputs. | [optional] |
| **workflowInputs5** | **String**| A fifth JSON or YAML file containing inputs. | [optional] |
| **workflowOptions** | **String**| JSON file containing configuration options for the execution of this workflow. | [optional] |
| **workflowType** | **String**| The workflow language for the file you submitted. Cromwell currently supports WDL and CWL. | [optional] [enum: WDL, CWL] |
| **workflowRoot** | **String**| The root object to be run. Only necessary for CWL submissions containing multiple objects (in an array). | [optional] |
| **workflowTypeVersion** | **String**| The specification version for the workflow language being used. For WDL, Cromwell currently supports draft-2 and 1.0. For CWL, Cromwell currently supports v1.0. | [optional] [enum: draft-2, 1.0, v1.0] |
| **labels** | **String**| JSON object of labels to apply to this workflow. | [optional] |
| **workflowDependencies** | **String**| ZIP file containing workflow source files that are used to resolve local imports. This zip bundle will be unpacked in a sandbox accessible to this workflow. | [optional] |
| **requestedWorkflowId** | **String**| An ID to ascribe to this workflow. Must be a JSON string in UUID-format. If not supplied a random ID will be generated for the workflow. | [optional] |

### Return type

[**WorkflowIdAndStatus**](WorkflowIdAndStatus.md)

### Authorization

[googleoauth](../README.md#googleoauth)

### HTTP request headers

 - **Content-Type**: multipart/form-data
 - **Accept**: application/json

### HTTP response details
| Status code | Description | Response headers |
|-------------|-------------|------------------|
| **201** | Successful Request |  -  |
| **400** | Invalid submission request |  -  |
| **500** | Internal Error |  -  |

<a name="submitBatch"></a>
# **submitBatch**
> List&lt;WorkflowIdAndStatus&gt; submitBatch(version, workflowInputs, workflowSource, workflowUrl, workflowOnHold, workflowOptions, workflowType, workflowTypeVersion, labels, workflowDependencies, requestedWorkflowId)

Submit a batch of workflows for execution

In instances where you want to run the same workflow multiple times with varying inputs you may submit a workflow batch. This endpoint is fundamentally the same as the standard submission endpoint with the exception that the inputs JSON will be an array of objects instead of a single object.

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.WorkflowsApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    WorkflowsApi apiInstance = new WorkflowsApi(defaultClient);
    String version = "v1"; // String | Cromwell API Version
    String workflowInputs = "workflowInputs_example"; // String | JSON file containing the inputs as an array of objects. Every element of the array will correspond to a single workflow. For WDL workflows a skeleton file can be generated from WOMtool using the \\\"inputs\\\" subcommand. When multiple files are specified, in case of key conflicts between multiple input JSON files, higher values of x in workflowInputs_x override lower values. For example, an input specified in workflowInputs_3 will override an input with the same name in workflowInputs or workflowInputs_2. Similarly, an input key specified in workflowInputs_5 will override an identical input key in any other input file.
    String workflowSource = "workflowSource_example"; // String | The workflow source file to submit for execution. Either workflow source or workflow url is required.
    String workflowUrl = "workflowUrl_example"; // String | URL which points to the workflow. Either workflow source or workflow url is required.
    Boolean workflowOnHold = true; // Boolean | Put workflow on hold upon submission. By default, it is taken as false.
    String workflowOptions = "workflowOptions_example"; // String | JSON file containing configuration options for the execution of this workflow.
    String workflowType = "WDL"; // String | The workflow language for the file you submitted. Cromwell currently supports WDL and CWL.
    String workflowTypeVersion = "draft-2"; // String | The specification version for the workflow language being used. For WDL, Cromwell currently supports draft-2 and 1.0. For CWL, Cromwell currently supports v1.0.
    String labels = "labels_example"; // String | JSON object of labels to apply to this workflow.
    String workflowDependencies = "workflowDependencies_example"; // String | ZIP file containing workflow source files that are used to resolve local imports. This zip bundle will be unpacked in a sandbox accessible to these workflows.
    String requestedWorkflowId = "requestedWorkflowId_example"; // String | A set of IDs to ascribe to these workflows. Must be a JSON list of strings in UUID-format. Must have the same number of entries and be in the same order as the workflow inputs list. If not supplied, random ID will be generated for the workflows.
    try {
      List<WorkflowIdAndStatus> result = apiInstance.submitBatch(version, workflowInputs, workflowSource, workflowUrl, workflowOnHold, workflowOptions, workflowType, workflowTypeVersion, labels, workflowDependencies, requestedWorkflowId);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling WorkflowsApi#submitBatch");
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
| **version** | **String**| Cromwell API Version | [default to v1] |
| **workflowInputs** | **String**| JSON file containing the inputs as an array of objects. Every element of the array will correspond to a single workflow. For WDL workflows a skeleton file can be generated from WOMtool using the \\\&quot;inputs\\\&quot; subcommand. When multiple files are specified, in case of key conflicts between multiple input JSON files, higher values of x in workflowInputs_x override lower values. For example, an input specified in workflowInputs_3 will override an input with the same name in workflowInputs or workflowInputs_2. Similarly, an input key specified in workflowInputs_5 will override an identical input key in any other input file. | |
| **workflowSource** | **String**| The workflow source file to submit for execution. Either workflow source or workflow url is required. | [optional] |
| **workflowUrl** | **String**| URL which points to the workflow. Either workflow source or workflow url is required. | [optional] |
| **workflowOnHold** | **Boolean**| Put workflow on hold upon submission. By default, it is taken as false. | [optional] |
| **workflowOptions** | **String**| JSON file containing configuration options for the execution of this workflow. | [optional] |
| **workflowType** | **String**| The workflow language for the file you submitted. Cromwell currently supports WDL and CWL. | [optional] [enum: WDL, CWL] |
| **workflowTypeVersion** | **String**| The specification version for the workflow language being used. For WDL, Cromwell currently supports draft-2 and 1.0. For CWL, Cromwell currently supports v1.0. | [optional] [enum: draft-2, 1.0, v1.0] |
| **labels** | **String**| JSON object of labels to apply to this workflow. | [optional] |
| **workflowDependencies** | **String**| ZIP file containing workflow source files that are used to resolve local imports. This zip bundle will be unpacked in a sandbox accessible to these workflows. | [optional] |
| **requestedWorkflowId** | **String**| A set of IDs to ascribe to these workflows. Must be a JSON list of strings in UUID-format. Must have the same number of entries and be in the same order as the workflow inputs list. If not supplied, random ID will be generated for the workflows. | [optional] |

### Return type

[**List&lt;WorkflowIdAndStatus&gt;**](WorkflowIdAndStatus.md)

### Authorization

[googleoauth](../README.md#googleoauth)

### HTTP request headers

 - **Content-Type**: multipart/form-data
 - **Accept**: application/json

### HTTP response details
| Status code | Description | Response headers |
|-------------|-------------|------------------|
| **200** | Successful Request |  -  |
| **400** | Malformed Workflow ID |  -  |
| **500** | Internal Error |  -  |

<a name="timing"></a>
# **timing**
> WorkflowIdAndStatus timing(version, id)

Get a visual diagram of a running workflow

Returns a javascript file which will render a Gantt chart for the requested workflow. The bars in the chart represent start and end times for individual task invocations. This javascript is intended to be embedded into another web page.

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.WorkflowsApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    WorkflowsApi apiInstance = new WorkflowsApi(defaultClient);
    String version = "v1"; // String | Cromwell API Version
    String id = "id_example"; // String | A workflow ID
    try {
      WorkflowIdAndStatus result = apiInstance.timing(version, id);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling WorkflowsApi#timing");
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
| **version** | **String**| Cromwell API Version | [default to v1] |
| **id** | **String**| A workflow ID | |

### Return type

[**WorkflowIdAndStatus**](WorkflowIdAndStatus.md)

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
| **404** | Workflow ID Not Found |  -  |
| **500** | Internal Error |  -  |

<a name="updateLabels"></a>
# **updateLabels**
> LabelsResponse updateLabels(version, id, labels)

Update labels for a workflow

Update multiple labels for an existing workflow. When supplying a label with a key unique to the workflow submission, a new label key/value entry is appended to that workflow&#39;s metadata. When supplying a label with a key that is already associated to the workflow submission, the original label value is updated with the new value for that workflow&#39;s metadata.

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.WorkflowsApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    WorkflowsApi apiInstance = new WorkflowsApi(defaultClient);
    String version = "v1"; // String | Cromwell API Version
    String id = "id_example"; // String | Workflow ID
    Object labels = null; // Object | Custom labels submitted as JSON. Example: {\"key-1\":\"value-1\",\"key-2\":\"value-2\"} 
    try {
      LabelsResponse result = apiInstance.updateLabels(version, id, labels);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling WorkflowsApi#updateLabels");
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
| **version** | **String**| Cromwell API Version | [default to v1] |
| **id** | **String**| Workflow ID | |
| **labels** | **Object**| Custom labels submitted as JSON. Example: {\&quot;key-1\&quot;:\&quot;value-1\&quot;,\&quot;key-2\&quot;:\&quot;value-2\&quot;}  | |

### Return type

[**LabelsResponse**](LabelsResponse.md)

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
| **403** | Workflow in terminal status |  -  |
| **404** | Workflow ID Not Found |  -  |
| **500** | Internal Error |  -  |

