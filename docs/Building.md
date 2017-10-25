_For the Doc-A-Thon_  
**Questions to answer and things to consider:**

1. Who is visiting the Building page?  
*Why would someone want to Build Cromwell?*
2. What do they need to know first?  
*someone is building Cromwell for the first time*
3. Is all the important information there? If not, add it!  
*add step by step instructions*
4. Are there things that don't need to be there? Remove them.  

5. Are the code and instructions accurate? Try it!

---
 **DELETE ABOVE ONCE COMPLETE**

---

`sbt assembly` will build a runnable JAR in `target/scala-2.12/`

Tests are run via `sbt test`.  Note that the tests do require Docker to be running.  To test this out while downloading the Ubuntu image that is required for tests, run `docker pull ubuntu:latest` prior to running `sbt test`