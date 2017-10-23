**Disclaimer**: This endpoint depends on hash values being published to the metadata, which only happens as of Cromwell 28.
Workflows run with prior versions of Cromwell cannot be used with this endpoint.
A `404 NotFound` will be returned when trying to use this endpoint if either workflow has been run on a prior version.

This endpoint returns the hash differences between 2 *completed* (successfully or not) calls.
The following query parameters are supported:

| Parameter |                                        Description                                        | Required |
|:---------:|:-----------------------------------------------------------------------------------------:|:--------:|
| workflowA | Workflow ID of the first call                                                             | yes      |
| callA     | Fully qualified name of the first call. **Including workflow name**. (see example below)  | yes      |
| indexA    | Shard index of the first call                                                             | depends  |
| workflowB | Workflow ID of the second call                                                            | yes      |
| callB     | Fully qualified name of the second call. **Including workflow name**. (see example below) | yes      |
| indexB    | Shard index of the second call                                                            | depends  |

About the `indexX` parameters: It is required if the call was in a scatter. Otherwise it should *not* be specified.
If an index parameter is wrongly specified, the call will not be found and the request will result in a 404 response.

cURL:

```
$ curl "http://localhost:8000/api/workflows/v1/callcaching/diff?workflowA=85174842-4a44-4355-a3a9-3a711ce556f1&callA=wf_hello.hello&workflowB=7479f8a8-efa4-46e4-af0d-802addc66e5d&callB=wf_hello.hello"
```

HTTPie:

```
$ http "http://localhost:8000/api/workflows/v1/callcaching/diff?workflowA=85174842-4a44-4355-a3a9-3a711ce556f1&callA=wf_hello.hello&workflowB=7479f8a8-efa4-46e4-af0d-802addc66e5d&callB=wf_hello.hello"
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 1274
Content-Type: application/json; charset=UTF-8
Date: Tue, 06 Jun 2017 16:44:33 GMT
Server: spray-can/1.3.3

{
  "callA": {
    "executionStatus": "Done",
    "workflowId": "85174842-4a44-4355-a3a9-3a711ce556f1",
    "callFqn": "wf_hello.hello",
    "jobIndex": -1,
    "allowResultReuse": true
  },
  "callB": {
    "executionStatus": "Done",
    "workflowId": "7479f8a8-efa4-46e4-af0d-802addc66e5d",
    "callFqn": "wf_hello.hello",
    "jobIndex": -1,
    "allowResultReuse": true
  },
  "hashDifferential": [
    {
      "hashKey": "command template",
      "callA": "4EAADE3CD5D558C5A6CFA4FD101A1486",
      "callB": "3C7A0CA3D7A863A486DBF3F7005D4C95"
    },
    {
      "hashKey": "input count",
      "callA": "C4CA4238A0B923820DCC509A6F75849B",
      "callB": "C81E728D9D4C2F636F067F89CC14862C"
    },
    {
      "hashKey": "input: String addressee",
      "callA": "D4CC65CB9B5F22D8A762532CED87FE8D",
      "callB": "7235E005510D99CB4D5988B21AC97B6D"
    },
    {
      "hashKey": "input: String addressee2",
      "callA": "116C7E36B4AE3EAFD07FA4C536CE092F",
      "callB": null
    }
  ]
}
```

The response is a JSON object with 3 fields:

- `callA` reports information about the first call, including its `allowResultReuse` value that will be used to determine whether or not this call can be cached to.
- `callB` reports information about the second call, including its `allowResultReuse` value that will be used to determine whether or not this call can be cached to.
- `hashDifferential` is an array in which each element represents a difference between the hashes of `callA` and `callB`.

*If this array is empty, `callA` and `callB` have the same hashes*.

Differences can be of 3 kinds:

- `callA` and `callB` both have the same hash key but their values are different.
For instance, in the example above, 

```json
{
  "hashKey": "input: String addressee",
  "callA": "D4CC65CB9B5F22D8A762532CED87FE8D",
  "callB": "7235E005510D99CB4D5988B21AC97B6D"
}
```

indicates that both `callA` and `callB` have a `String` input called `addressee`, but different values were used at runtime, resulting in different MD5 hashes.

- `callA` has a hash key that `callB` doesn't have
For instance, in the example above, 

```json
{
  "hashKey": "input: String addressee2",
  "callA": "116C7E36B4AE3EAFD07FA4C536CE092F",
  "callB": null
}
```

indicates that `callA` has a `String` input called `addressee2` that doesn't exist in `callB`. For that reason the value of the second field is `null`.

- `callB` has a hash key that `callA` doesn't have. This is the same case as above but reversed.

If no cache entry for `callA` or `callB` can be found, the response will be in the following format:

```
HTTP/1.1 404 NotFound
Content-Length: 178
Content-Type: application/json; charset=UTF-8
Date: Tue, 06 Jun 2017 17:02:15 GMT
Server: spray-can/1.3.3

{
  "status": "error",
  "message": "Cannot find call 479f8a8-efa4-46e4-af0d-802addc66e5d:wf_hello.hello:-1"
}
```

If neither `callA` nor `callB` can be found, the response will be in the following format:


```
HTTP/1.1 404 NotFound
Content-Length: 178
Content-Type: application/json; charset=UTF-8
Date: Tue, 06 Jun 2017 17:02:15 GMT
Server: spray-can/1.3.3

{
  "status": "error",
  "message": "Cannot find calls 5174842-4a44-4355-a3a9-3a711ce556f1:wf_hello.hello:-1, 479f8a8-efa4-46e4-af0d-802addc66e5d:wf_hello.hello:-1"
}
```

If the query is malformed and required parameters are missing, the response will be in the following format:

```
HTTP/1.1 400 BadRequest
Content-Length: 178
Content-Type: application/json; charset=UTF-8
Date: Tue, 06 Jun 2017 17:02:15 GMT
Server: spray-can/1.3.3
{
  "status": "fail",
  "message": "Wrong parameters for call cache diff query:\nmissing workflowA query parameter\nmissing callB query parameter",
  "errors": [
    "missing workflowA query parameter",
    "missing callB query parameter"
  ]
}
```