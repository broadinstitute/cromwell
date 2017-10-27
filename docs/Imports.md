Cromwell implements imports in both [run](CommandLine) and [server](api/POST_api_workflows_version) mode.
You can either include imported files when submitting the workflow, or reference them via `http` or `https`.

To include your own files, supply a .zip file that contains all required dependencies.
The directory structure of the .zip file should match the import paths in the workflow source.

Here's an example workflow in WDL:

_workflow.wdl_
```wdl
import "https://github.com/broadinstitute/cromwell/blob/master/engine/src/main/resources/3step.wdl" as http_import
import "imports/imported.wdl" as provided_import

workflow my_workflow {
    ...
}
```

_imports.zip_
```
imports
└── imported.wdl
```

---

In Run mode, a sample command to run _workflow.wdl_ would be  
```$ java -jar cromwell.jar run workflow.wdl --imports imports.zip```


In Server Mode, pass in a .zip file using the parameter `workflowDependencies` via the [submit](api/POST_api_workflows_version) endpoint.
