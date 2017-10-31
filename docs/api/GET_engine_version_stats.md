This endpoint returns some basic statistics on the current state of the engine. At the moment that includes the number of running workflows and the number of active jobs. 

cURL:
```
$ curl http://localhost:8000/engine/v1/stats
```

HTTPie:
```
$ http http://localhost:8000/engine/v1/stats
```

Response:
```
"date": "Sun, 18 Sep 2016 14:38:11 GMT",
"server": "spray-can/1.3.3",
"content-length": "33",
"content-type": "application/json; charset=UTF-8"

{
  "workflows": 3,
  "jobs": 10
}
```