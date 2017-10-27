_For the Doc-A-Thon_  
**Questions to answer and things to consider:**

1. Who is visiting the Imports page?  
*did they come from the Command line page? from the API docs?*
2. What do they need to know first?  
*What are imports? Why might they want to use them? Common use cases?*
3. Is all the important information there? If not, add it!  

4. Are there things that don't need to be there? Remove them.  
*What is Single Workflow Runner? It doesn't exist anywhere else in the docs. Either explain it or remove it!*
5. Are the code and instructions accurate? Try it!

---
 **DELETE ABOVE ONCE COMPLETE**

---


Import statements inside of a workflow file are supported by Cromwell when running in Server mode as well as Single Workflow Runner Mode.

In Single Workflow Runner Mode, you pass in a zip file which includes the WDL files referenced by the import statements. Cromwell requires the zip file to be passed in as a command line argument, as explained by the section [Run](CommandLine).

For example, given a workflow `wf.wdl` and an imports directory `WdlImports.zip`, a sample command would be:
```
java -jar cromwell.jar wf.wdl wf.inputs - - WdlImports.zip
```

In Server Mode, you pass in a zip file using the parameter `workflowDependencies` via the [POST /api/workflows/:version](/restapi#post-apiworkflowsversion) endpoint.

**Imports via HTTPS**

# TBC by Kcibul