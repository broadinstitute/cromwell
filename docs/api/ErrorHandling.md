Requests that Cromwell can't process return a failure in the form of a JSON response respecting the following JSON schema:

```
{
  "$schema": "http://json-schema.org/draft-04/schema#",
  "description": "Error response schema",
  "type": "object",
  "properties": {
    "status": {
      "enum": [ "fail", "error"]
    },
    "message": {
      "type": "string"
    },
    "errors": {
      "type": "array",
      "minItems": 1,
      "items": { "type": "string" },
      "uniqueItems": true
    }
  },
  "required": ["status", "message"]
}
```

The `status` field can be `"fail"` or `"error"`.  

`"fail"` means that the request was invalid and/or data validation failed. `"fail"` status is most likely returned with a 4xx HTTP Status code.  

*For example,*

```
{
  "status": "fail",
  "message": "Workflow input processing failed.",
  "errors": [
    "Required workflow input 'helloworld.input' not specified."
  ]
}
```

`"error"` means that an error occurred while processing the request. `"error"` status is most likely returned with a 5xx HTTP Status code.  

*For example,*

```
{
  "status": "error",
  "message": "Connection to the database failed."
}
```

The `message` field contains a short description of the error.

The `errors` field is optional and may contain additional information about why the request failed.