# Workflow Options Overview

Workflow options can affect the execution of a single workflow without having to change configuration options or restart Cromwell. 

Unless otherwise specified you can expect workflow options to override any hard-coded defaults in Cromwell or defaults provided in the [configuration file](../Configuring), but to be overridden by any values provided in the workflow definition file itself (WDL or CWL).

Provide workflow options via a JSON file that toggles various options for running the workflow. This can be supplied at workflow-submit time either via the (**TODO: make these links**) CLI or the REST endpoint.

Example:

```json
{
  "jes_gcs_root": "gs://my-bucket/workflows",
  "google_project": "my_google_project",
  "refresh_token": "1/Fjf8gfJr5fdfNf9dk26fdn23FDm4x"
}
```

Several workflow options apply to every task, regardless of backend: [Global Options](Global).

Some workflow options apply only to tasks running on the Google Pipelines API backend: [Google Options](Google)
