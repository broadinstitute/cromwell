This endpoint is used to update multiple labels for an existing workflow. When supplying a label with a key unique to the workflow submission, a new label key/value entry is appended to that workflow's metadata. When supplying a label with a key that is already associated to the workflow submission, the original label value is updated with the new value for that workflow's metadata.

The [labels](/labels) must be a mapping of key/value pairs in JSON format that are sent via the PATCH body. The request content type must be
`application/json`.

cURL:

```
$ curl -X PATCH --header "Content-Type: application/json" -d "{\"label-key-1\":\"label-value-1\", \"label-key-2\": \"label-value-2\"}" "http://localhost:8000/api/workflows/v1/c4c6339c-8cc9-47fb-acc5-b5cb8d2809f5/labels"
```

HTTPie:

```
$ echo '{"label-key-1":"label-value-1", "label-key-2": "label-value-2"}' | http PATCH "http://localhost:8000/api/workflows/v1/c4c6339c-8cc9-47fb-acc5-b5cb8d2809f5/labels"
```

Response:
```
{ "id": "c4c6339c-8cc9-47fb-acc5-b5cb8d2809f5",
  "labels":
    {
      "label-key-1": "label-value-1",
      "label-key-2": "label-value-2"
    }
}
```