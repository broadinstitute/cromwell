_For the Doc-A-Thon_  
**Questions to answer and things to consider:**

1. Who is visiting the API/backends page?  

2. What do they need to know first?  

3. Is all the important information there? If not, add it!  

4. Are there things that don't need to be there? Remove them.  

5. Are the code and instructions accurate? Try it!

---
 **DELETE ABOVE ONCE COMPLETE**

---


This endpoint returns a list of the backends supported by the server as well as the default backend.

cURL:
```
$ curl http://localhost:8000/api/workflows/v1/backends
```

HTTPie:
```
$ http http://localhost:8000/api/workflows/v1/backends
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 379
Content-Type: application/json; charset=UTF-8
Date: Mon, 03 Aug 2015 17:11:28 GMT
Server: spray-can/1.3.3

{
  "supportedBackends": ["JES", "LSF", "Local", "SGE"],
  "defaultBackend": "Local"
}
```