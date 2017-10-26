_For the Doc-A-Thon_  
**Questions to answer and things to consider:**

1. Who is visiting the API/outputs page?  

2. What do they need to know first?  

3. Is all the important information there? If not, add it!  
*Add a sentence or two at the beginning to explain the output (what they'll see, why it's useful, why they care)*
4. Are there things that don't need to be there? Remove them.  

5. Are the code and instructions accurate? Try it!

---
 **DELETE ABOVE ONCE COMPLETE**

---

cURL:

```
$ curl http://localhost:8000/api/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/outputs
```

HTTPie:

```
$ http http://localhost:8000/api/workflows/v1/e442e52a-9de1-47f0-8b4f-e6e565008cf1/outputs
```

Response:
```
HTTP/1.1 200 OK
Content-Length: 241
Content-Type: application/json; charset=UTF-8
Date: Thu, 04 Jun 2015 12:15:33 GMT
Server: spray-can/1.3.3

{
    "id": "e442e52a-9de1-47f0-8b4f-e6e565008cf1",
    "outputs": {
        "three_step.cgrep.count": 8,
        "three_step.ps.procs": "/var/folders/kg/c7vgxnn902lc3qvc2z2g81s89xhzdz/T/stdout2814345504446060277.tmp",
        "three_step.wc.count": 8
    }
}
```