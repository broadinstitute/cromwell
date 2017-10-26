_For the Doc-A-Thon_  
**Questions to answer and things to consider:**

1. Who is visiting the API/Stats page?  

2. What do they need to know first?  

3. Is all the important information there? If not, add it!  
*Why would they need these stats?*
4. Are there things that don't need to be there? Remove them.  

5. Are the code and instructions accurate? Try it!

---
 **DELETE ABOVE ONCE COMPLETE**

---

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