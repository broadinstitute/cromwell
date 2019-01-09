Sometimes you might want to break up a long WDL file into smaller components for easier maintenance. Sometimes you'll want to do it to reuse components in multiple workflows.  Have no fear, imports are here!  Imports allow you to reference other WDL files that contain entire workflows or even just raw tasks.

To import a WDL, you can use the `import` WDL construct at the top of your workflow:

```
import "<resource>" as <alias>
```

There are two types of resources that are supported in imports: *http(s)* and *file-path based*.  Any public http(s) based URL can be used as the resource for an import, such as a website, github or a GA4GH compliant TES endpoint.  For example:

```wdl
import "http://mywdlrepository/my.wdl" as http_import1
import "https://github.com/broadinstitute/cromwell/blob/master/engine/src/main/resources/3step.wdl" as http_import2
```
To use a file-based import resource, provide a ZIP bundle of your resources and then use a path relative to that ZIP in your import statement. For example:

```wdl
import "my-wdl-in-the-root-directory.wdl" as file_import1
import "my/wdl/sub/directory/example.wdl" as file_import2
```

Imports from your submitted workflow are evaluated relative to the base of the zip file. In other cases, import paths are relative to the file you are currently importing from, for example:

>If there exists a `my/wdl/sub/directory/imports/importing_an_import.wdl`, `my/wdl/sub/directory/example.wdl` could import it relatively like so:
>```wdl
>import "imports/importing_an_import.wdl" as file_import3
>```

Here's a complete example showing both http(s) and file-based imports workflow in WDL:

_workflow.wdl_
```wdl
import "https://github.com/broadinstitute/cromwell/blob/master/engine/src/main/resources/3step.wdl" as http_import
import "imports/imported.wdl" as provided_import

workflow my_workflow {
    ...
}
```

Note that we said "`import "imports/imported.wdl`" in the workflow so we must have a `imports/imported.wdl` structure in the imports file:
_imports.zip_
```
imports
└── imported.wdl
```

A more common scenario might have the imports at the root of the imports.zip:

_workflow.wdl_
```wdl
import "my_wdl_1.wdl"
import "my_wdl_2.wdl"
```
_imports.zip_
```
my_wdl_1.wdl
my_wdl_2.wdl
```

---
The mechanism to provide the ZIP file of resources to be imported differs between [Run](Modes#run) and [Server](Modes#server) mode.

In [Run](Modes#run) mode, a sample command to run _workflow.wdl_ would be:  
```
$ java -jar cromwell.jar run workflow.wdl --imports imports.zip
```

In [Server](Modes#server) mode, pass in a ZIP file using the parameter `workflowDependencies` via the [Submit](api/RESTAPI#submit-a-workflow-for-execution) endpoint.

