_For the Doc-A-Thon_  
**Questions to answer and things to consider:**

1. Who is visiting the API/Query page?  
*Why aren't they visiting the GET API Query page? Maybe they came there from there* 
2. What do they need to know first?  

3. Is all the important information there? If not, add it!  

4. Are there things that don't need to be there? Remove them.  

5. Are the code and instructions accurate? Try it!

---
 **DELETE ABOVE ONCE COMPLETE**

---

This endpoint allows for querying workflows based on the same criteria
as [GET /api/workflows/:version/query](/restpi#get-apiworkflowsversionquery).

Instead of specifying query parameters in the URL, the parameters
must be sent via the POST body. The request content type must be
`application/json`. The json should be a list of objects. Each json
object should contain a different criterion.

cURL:

```
$ curl -X POST --header "Content-Type: application/json" -d "[{\"start\": \"2015-11-01T00:00:00-04:00\"}, {\"end\": \"2015-11-04T00:00:00-04:00\"}, {\"status\": \"Failed\"}, {\"status\": \"Succeeded\"}, {\"page\": \"1\"}, {\"pagesize\": \"10\"}]" "http://localhost:8000/api/workflows/v1/query"
```

HTTPie:

```
$ echo "[{\"start\": \"2015-11-01T00:00:00-04:00\"}, {\"end\": \"2015-11-04T00:00:00-04:00\"}, {\"status\": \"Failed\"}, {\"status\": \"Succeeded\"}, {\"page\": \"1\"}, {\"pagesize\": \"10\"}]" | http "http://localhost:8000/api/workflows/v1/query"
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 133
Content-Type: application/json; charset=UTF-8
Date: Tue, 02 Jun 2015 18:06:56 GMT
Server: spray-can/1.3.3

{
  "results": [
    {
      "name": "w",
      "id": "fdfa8482-e870-4528-b639-73514b0469b2",
      "status": "Succeeded",
      "end": "2015-11-01T07:45:52.000-05:00",
      "start": "2015-11-01T07:38:57.000-05:00"
    },
    {
      "name": "hello",
      "id": "e69895b1-42ed-40e1-b42d-888532c49a0f",
      "status": "Succeeded",
      "end": "2015-11-01T07:45:30.000-05:00",
      "start": "2015-11-01T07:38:58.000-05:00"
    },
    {
      "name": "crasher",
      "id": "ed44cce4-d21b-4c42-b76d-9d145e4d3607",
      "status": "Failed",
      "end": "2015-11-01T07:45:44.000-05:00",
      "start": "2015-11-01T07:38:59.000-05:00"
    }
  ]
}
```