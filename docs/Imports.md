Import statements inside of a workflow file are supported by Cromwell when running in Server mode as well as Single Workflow Runner Mode.

In Single Workflow Runner Mode, you pass in a zip file which includes the WDL files referenced by the import statements. Cromwell requires the zip file to be passed in as a command line argument, as explained by the section [run](#run).

For example, given a workflow `wf.wdl` and an imports directory `WdlImports.zip`, a sample command would be:
```
java -jar cromwell.jar wf.wdl wf.inputs - - WdlImports.zip
```

In Server Mode, you pass in a zip file using the parameter `workflowDependencies` via the [POST /api/workflows/:version](#post-apiworkflowsversion) endpoint.