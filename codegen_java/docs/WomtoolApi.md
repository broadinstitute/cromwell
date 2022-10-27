# WomtoolApi

All URIs are relative to *http://localhost*

| Method | HTTP request | Description |
|------------- | ------------- | -------------|
| [**describe**](WomtoolApi.md#describe) | **POST** /api/womtool/{version}/describe | Machine-readable description of a workflow, including inputs and outputs |


<a name="describe"></a>
# **describe**
> WorkflowDescription describe(version, workflowSource, workflowUrl, workflowInputs, workflowType, workflowTypeVersion)

Machine-readable description of a workflow, including inputs and outputs

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.WomtoolApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    WomtoolApi apiInstance = new WomtoolApi(defaultClient);
    String version = "v1"; // String | Cromwell API Version
    String workflowSource = "workflowSource_example"; // String | The workflow source file to submit for execution. Either workflow source or workflow url is required.
    String workflowUrl = "workflowUrl_example"; // String | URL which points to the workflow. Either workflow source or workflow url is required.
    String workflowInputs = "workflowInputs_example"; // String | JSON or YAML file containing the inputs as an object.
    String workflowType = "WDL"; // String | The workflow language for the file you submitted. Cromwell currently supports WDL and CWL.
    String workflowTypeVersion = "draft-2"; // String | The specification version for the workflow language being used. For WDL, Cromwell currently supports draft-2 and 1.0. For CWL, Cromwell currently supports v1.0.
    try {
      WorkflowDescription result = apiInstance.describe(version, workflowSource, workflowUrl, workflowInputs, workflowType, workflowTypeVersion);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling WomtoolApi#describe");
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
| **workflowInputs** | **String**| JSON or YAML file containing the inputs as an object. | [optional] |
| **workflowType** | **String**| The workflow language for the file you submitted. Cromwell currently supports WDL and CWL. | [optional] [enum: WDL, CWL] |
| **workflowTypeVersion** | **String**| The specification version for the workflow language being used. For WDL, Cromwell currently supports draft-2 and 1.0. For CWL, Cromwell currently supports v1.0. | [optional] [enum: draft-2, 1.0, v1.0] |

### Return type

[**WorkflowDescription**](WorkflowDescription.md)

### Authorization

[googleoauth](../README.md#googleoauth)

### HTTP request headers

 - **Content-Type**: multipart/form-data
 - **Accept**: application/json

### HTTP response details
| Status code | Description | Response headers |
|-------------|-------------|------------------|
| **200** | Workflow description. |  -  |

