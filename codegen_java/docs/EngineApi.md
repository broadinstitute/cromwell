# EngineApi

All URIs are relative to *http://localhost*

| Method | HTTP request | Description |
|------------- | ------------- | -------------|
| [**engineStatus**](EngineApi.md#engineStatus) | **GET** /engine/{version}/status | Return the current health status of any monitored subsystems |
| [**engineVersion**](EngineApi.md#engineVersion) | **GET** /engine/{version}/version | Return the version of this Cromwell server |


<a name="engineStatus"></a>
# **engineStatus**
> StatusResponse engineStatus(version)

Return the current health status of any monitored subsystems

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.EngineApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    EngineApi apiInstance = new EngineApi(defaultClient);
    String version = "v1"; // String | Cromwell API Version
    try {
      StatusResponse result = apiInstance.engineStatus(version);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling EngineApi#engineStatus");
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

[**StatusResponse**](StatusResponse.md)

### Authorization

[googleoauth](../README.md#googleoauth)

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

### HTTP response details
| Status code | Description | Response headers |
|-------------|-------------|------------------|
| **200** | All subsystems report an \&quot;ok\&quot; status |  -  |
| **500** | At least one subsystem does not have an \&quot;ok\&quot; status |  -  |

<a name="engineVersion"></a>
# **engineVersion**
> VersionResponse engineVersion(version)

Return the version of this Cromwell server

### Example
```java
// Import classes:
import cromwell.client.ApiClient;
import cromwell.client.ApiException;
import cromwell.client.Configuration;
import cromwell.client.auth.*;
import cromwell.client.models.*;
import cromwell.client.api.EngineApi;

public class Example {
  public static void main(String[] args) {
    ApiClient defaultClient = Configuration.getDefaultApiClient();
    defaultClient.setBasePath("http://localhost");
    
    // Configure OAuth2 access token for authorization: googleoauth
    OAuth googleoauth = (OAuth) defaultClient.getAuthentication("googleoauth");
    googleoauth.setAccessToken("YOUR ACCESS TOKEN");

    EngineApi apiInstance = new EngineApi(defaultClient);
    String version = "v1"; // String | Cromwell API Version
    try {
      VersionResponse result = apiInstance.engineVersion(version);
      System.out.println(result);
    } catch (ApiException e) {
      System.err.println("Exception when calling EngineApi#engineVersion");
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

[**VersionResponse**](VersionResponse.md)

### Authorization

[googleoauth](../README.md#googleoauth)

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

### HTTP response details
| Status code | Description | Response headers |
|-------------|-------------|------------------|
| **200** | Successful Request |  -  |

