_For the Doc-A-Thon_  
**Questions to answer and things to consider:**

1. Who is visiting the API/Cromwell version page?  

2. What do they need to know first?  
*Right now we just have one version, right? Is this a useful endpoint?*
3. Is all the important information there? If not, add it!  

4. Are there things that don't need to be there? Remove them.  

5. Are the code and instructions accurate? Try it!

---
 **DELETE ABOVE ONCE COMPLETE**

---

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