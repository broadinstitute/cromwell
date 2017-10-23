This endpoint returns the version of the Cromwell engine.

cURL:
```
$ curl http://localhost:8000/engine/v1/version
```

HTTPie:
```
$ http http://localhost:8000/engine/v1/version
```

Response:
```
"date": "Sun, 18 Sep 2016 14:38:11 GMT",
"server": "spray-can/1.3.3",
"content-length": "33",
"content-type": "application/json; charset=UTF-8"

{
  "cromwell": 23-8be799a-SNAP
}
```