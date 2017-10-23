cURL:

```
$ curl http://localhost:8000/api/workflows/v1/69d1d92f-3895-4a7b-880a-82535e9a096e/status
```

HTTPie:

```
$ http http://localhost:8000/api/workflows/v1/69d1d92f-3895-4a7b-880a-82535e9a096e/status
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 74
Content-Type: application/json; charset=UTF-8
Date: Tue, 02 Jun 2015 18:06:56 GMT
Server: spray-can/1.3.3

{
    "id": "69d1d92f-3895-4a7b-880a-82535e9a096e",
    "status": "Succeeded"
}
```