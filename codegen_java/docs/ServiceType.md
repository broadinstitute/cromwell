

# ServiceType

Type of a GA4GH service

## Properties

| Name | Type | Description | Notes |
|------------ | ------------- | ------------- | -------------|
|**group** | **String** | Namespace in reverse domain name format. Use &#x60;org.ga4gh&#x60; for implementations compliant with official GA4GH specifications. For services with custom APIs not standardized by GA4GH, or implementations diverging from official GA4GH specifications, use a different namespace (e.g. your organization&#39;s reverse domain name). |  |
|**artifact** | **String** | Name of the API or GA4GH specification implemented. Official GA4GH types should be assigned as part of standards approval process. Custom artifacts are supported. |  |
|**version** | **String** | Version of the API or specification. GA4GH specifications use semantic versioning. |  |



