Sometimes you might want to break up your 1000 line WDL file into smaller components for easier maintenance or for reuse in multiple workflows.  Have no fear, imports are here to help!  Imports allow you to reference other WDL files that contain entire workflows or even just raw tasks.

To import a WDL, you can use the `import` WDL construct at the top of your workflow

```
import "<resource>" as <alias>
```

There are two types of resources that are supported in imports: *http(s)* and *file-path based*.  Any public http(s) based URL can be used as the resource for an import, such as a website, github or a GA4GH compliant TES endpoint.  For example:

```wdl
import "http://mywdlrepository/my.wdl" as http_import1
import "https://github.com/broadinstitute/cromwell/blob/master/engine/src/main/resources/3step.wdl" as http_import2
```
To use a file-based import resource, you must provide a ZIP bundle of your resources and then use a path relative to that ZIP in your import statement. For example:

```wdl
import "my-wdl-in-the-root-directory.wdl" as file_import1
import "my/wdl/sub/directory/example.wdl" as file_import2
```

Here's a complete example showing both http(s) and file-based imports workflow in WDL:

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
The mechanism to provide the ZIP file of resources to be imported differ between Run and Server mode.

In Run mode, a sample command to run _workflow.wdl_ would be  
```$ java -jar cromwell.jar run workflow.wdl --imports imports.zip```

In Server Mode, pass in a .zip file using the parameter `workflowDependencies` via the [submit](api/POST_api_workflows_version) endpoint.

