
For built-in documentation of Cromwell command line usage, run the Cromwell JAR file with no arguments:

`$ java -jar cromwell-<versionNumber>.jar`

You will get a usage message like the following:

```bash
cromwell 29
Usage: java -jar /path/to/cromwell.jar [server|run] [options] <args>...

  --help                   Cromwell - Workflow Execution Engine
  --version                
Command: server
Starts a web server on port 8000.  See the web server documentation for more details about the API endpoints.
Command: run [options] workflow-source
Run the workflow and print out the outputs in JSON format.
  workflow-source          Workflow source file or workflow url .
  --workflow-root <value>  Workflow root
  -i, --inputs <value>     Workflow inputs file.
  -o, --options <value>    Workflow options file.
  -t, --type <value>       Workflow type.
  -v, --type-version <value>
                           Workflow type version.
  -l, --labels <value>     Workflow labels file.
  -p, --imports <value>    A directory or zipfile to search for workflow imports.
  -m, --metadata-output <value>
                           An optional JSON file path to output metadata.
```

Cromwell's Server and Run modes can be invoked with the `server` and `run` arguments respectively. More information on these Cromwell modes can be found in [Modes](Modes).

The Cromwell jar file can be built as described in [Building](Building). 

## `server`

`server` mode accepts no arguments and runs Cromwell as a web server that accepts REST requests. See the documentation for [Cromwell's REST endpoints](/api/RESTAPI) for how to interact with Cromwell in `server` mode.

## `run`

`run` mode executes a single workflow in Cromwell and then exits.

* **`workflow-source`**  
The single required argument. It can be either a local path or a remote URL pointing to the workflow source file.
 
* **`--inputs`**  
An optional file of workflow inputs.  Although optional, it is a best practice to use an inputs file to satisfy workflow
requirements rather than hardcoding inputs directly into a workflow source file.

* **`--options`**  
An optional file of workflow options.  Some options are global (supported by all backends), while others are backend-specific. See the [Workflow Options](wf_options/Overview) for more details.

* **`--type`**  
An optional parameter to specify the language for the workflow source. As of Cromwell 29 any value specified for this parameter is currently ignored and internally the value `WDL` is used.

* **`--type-version`**  
An optional parameter to specify the version of the language for the workflow source. Currently any specified value is ignored.

* **`--labels`**  
An optional parameter to specify a file of JSON key-value label pairs to associate with the workflow.

* **`--imports`**  
You have the option of importing WDL workflows or tasks to use within your workflow, known as sub-workflows.
If you use sub-workflows within your primary workflow then you must include a ZIP file with the WDL import files.
See the documentation on [Imports](Imports) for more information.

* **`--metadata-output`**  
You can specify a filename where Cromwell will write workflow metadata JSON such as start/end timestamps, status, inputs and outputs. By default Cromwell does not write workflow metadata. The metadata format in the `--metadata-output` file is the same as described for the [REST API](api/RESTAPI#get-workflow-and-call-level-metadata-for-a-specified-workflow).

* **`--version`**  
The `--version` option prints the version of Cromwell and exits.

* **`--help`**  
The `--help` option prints the full help text above and exits.
