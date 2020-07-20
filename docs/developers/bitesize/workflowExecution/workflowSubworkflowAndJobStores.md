# Workflow Execution: Database Tables: Workflow, Subworkflow, and Job Stores

## Database Tables: Workflow, Subworkflow, and Job Stores

Cromwell uses the workflow, subworkflow and job store tables to hold data related to submitted or running workflows.
Once a workflow reaches a terminal state all data for that workflow should be deleted from these tables.

### Workflow Store / `WORKFLOW_STORE_ENTRY`

`WORKFLOW_STORE_ENTRY` holds data received in a workflow submission (workflow sources, inputs, options etc.)
and workflow-scoped execution data (e.g. submission time, status, fields to support
running [Horizontal Cromwell](../horicromtal.md) etc).

### Job Store / `JOB_STORE_ENTRY`

`JOB_STORE_ENTRY` holds data for *completed* jobs within a workflow. Jobs that are still running or have not yet been
started will not have rows in this table. The main purpose of the job store table is to support resuming execution of
a workflow when Cromwell is restarted by recovering the outputs of completed jobs. This table is closely related to
`JOB_STORE_SIMPLETON_ENTRY` which holds the [simpleton](../general/simpletons.md) values comprising a job's outputs,
and loosely related to the [job key/value store (`JOB_KEY_VALUE_ENTRY`)](jobKeyValueStore.md) which holds other
job-scoped data important in recovering jobs on Cromwell restart.

### Subworkflow Store / `SUB_WORKFLOW_STORE_ENTRY`

`SUB_WORKFLOW_STORE_ENTRY` holds data for subworkflows that have begun execution. The rows in this table persist the fact
that particular subworkflows corresponding to a call FQN and index were started and assigned a workflow ID.
The completed jobs within these subworkflows will be recorded in the job store described above, linking to the
subworkflows in this table by the subworkflow's ID.
