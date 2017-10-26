_For the Doc-A-Thon_  
**Questions to answer and things to consider:**

1. Who is visiting the API/engine status page?  

2. What do they need to know first?  

3. Is all the important information there? If not, add it!  
*When would someone want this information?*
4. Are there things that don't need to be there? Remove them.  

5. Are the code and instructions accurate? Try it!

---
 **DELETE ABOVE ONCE COMPLETE**

---

Cromwell will track the status of underlying systems on which it depends, typically connectivity to the database and to Docker Hub. The `status` endpoint will return the current status of these systems. 

Response:

This endpoint will return an `Internal Server Error` if any systems are marked as failing or `OK` otherwise. The exact response will vary based on what is being monitored but adheres to the following pattern. Each system has a boolean `ok` field and an optional array field named `messages` which will be populated with any known errors if the `ok` status is `false`:

```
{
  "System Name 1": {
    "ok": false,
    "messages": [
      "Unknown status"
    ]
  },
  "System Name 2": {
    "ok": true
  }
}

```