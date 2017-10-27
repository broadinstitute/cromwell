Cromwell implements imports in both [run](CommandLine) and [server](api/POST_api_workflows_version) mode.
Imported files can be provided when submitting the workflow or referenced via `http` or `https`.
In order to provide your own files, Cromwell lets you supply a zip file that should contain all required dependencies.
The directory structure of the zip file should match the import paths in the workflow source.

For example, in WDL:

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
```java -jar cromwell.jar run workflow.wdl --imports imports.zip```


In Server Mode, pass in a zip file using the parameter `workflowDependencies` via the [submit](/restapi#post-apiworkflowsversion) endpoint.
