

# Service

GA4GH service

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



