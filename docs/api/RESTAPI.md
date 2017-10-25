_For the Doc-A-Thon_  
**Questions to answer and things to consider:**

1. Who is visiting the API page?  

2. What do they need to know first?  

3. Is all the important information there? If not, add it!  

4. Are there things that don't need to be there? Remove them.  
*Does anyone use HTTPie? Jeff says no, if you disagree quabble with him*
5. Are the code and instructions accurate? Try it!

---
 **DELETE ABOVE ONCE COMPLETE**

---


The `server` subcommand on the executable JAR will start an HTTP server which can accept workflow files to run as well as check status and output of existing workflows.

The following sub-sections define which HTTP Requests the web server can accept and what they will return. Example HTTP requests are given in [HTTPie](https://github.com/jkbrzt/httpie) and [cURL](https://curl.haxx.se/)

**REST API Versions**

All web server requests include an API version in the url. The current version is `v1`.