# Workflow Execution: Major Actors

* **Word Count:** 245

## Major Actor Hierarchy

At the highest level, these are the main actors involved in workflow execution.

![high level overview diagram](WorkflowExecutionHighLevelOverview.png)

## Actors and their Purposes

### WorkflowManagerActor

The `WorkflowManagerActor` is responsible for:

* Polling the `WorkflowStore` at pre-configured intervals.
* Starting new workflows
* Tracking, supervising and aborting running workflows
* Parent actor for all `WorkflowActor`s

### WorkflowActor(s)

The `WorkflowActor` is responsible for:
 
* Co-ordinating the stages of a workflow lifecycle from parsing through to finalization.
* Parent actor of the `WorkflowExecutionActor` which runs the workflow's jobs.

### WorkflowExecutionActor(s)

The `WorkflowExecutionActor` is responsible for:

* Starting jobs and sub-workflows as soon as they are able to run.
    * Based on values in the (in-memory) ValueStore and ExecutionStore objects.
* Parent actor for all `EngineJobExecutionActor`s and `SubWorkflowExecutionActor`s.

### EngineJobExecutionActor(s)

Each `EngineJobExecutionActor` (EJEA) is responsible for:

* Running a single job.
    * A "job" is a command line instruction to run on a backend.
    * Multiple shards for a single call each get their own EJEA.
    * Multiple attempts to run the same job operate within the same EJEA
* Respects hog-limiting
* Checks the call cache and job store to avoid running the job if it doesn't have to.
* Triggers job initialization, execution and finalization at appropriate times.

### SubWorkflowExecutionActor(s)

Each `SubWorkflowExecutionActor` is responsible for:

* Running a single sub-workflow.
* Parent actor for a new `WorkflowExecutionActor` (see above) created to run the sub-workflow.

## Major Actor Hierarchy (in context)

The above diagram omitted a lot of details. This diagram attempts to show a little more of the
context:

![high level overview in context diagram](WorkflowExecutionHighLevelOverviewInContext.png)

## See Also 

* EngineJobExecutionActor (**TODO**)
* Backend Execution Actors (**TODO**)
