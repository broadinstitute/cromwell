_For the Doc-A-Thon_  
**Questions to answer and things to consider:**

1. Who is visiting the Errors page?  

2. What do they need to know first?  
*What is the purpose of this page?*
3. Is all the important information there? If not, add it!  
*Frequently seen errors?*
4. Are there things that don't need to be there? Remove them.  

5. Are the code and instructions accurate? Try it!

---
 **DELETE ABOVE ONCE COMPLETE**

---


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

The `status` field can take two values:
> "fail" means that the request was invalid and/or data validation failed. "fail" status is most likely returned with a 4xx HTTP Status code.
e.g.

```
{
  "status": "fail",
  "message": "Workflow input processing failed.",
  "errors": [
    "Required workflow input 'helloworld.input' not specified."
  ]
}
```

> "error" means that an error occurred while processing the request. "error" status is most likely returned with a 5xx HTTP Status code.
e.g.

```
{
  "status": "error",
  "message": "Connection to the database failed."
}
```

The `message` field contains a short description of the error.

The `errors` field is optional and may contain additional information about why the request failed.