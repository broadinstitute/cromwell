# Workflow Execution: Execution Store and Value Store Examples

## Introduction

This page provides run-throughs to give insight into how the
[Execution Store](executionStore.md) and [Value Store](valueStore.md) work in
practice in some example situations. 

## Handling a single task call

To begin, consider this simple workflow. It has a single task call whose
result is exposed as an output String:

```wdl
version 1.0

workflow single_task_workflow {
  call single_task
  
  output {
    String string_out = single_task.string_out
  }
}

task single_task {
  command {
    echo hello
  }
  output {
    String string_out = "hello"
  }
}
```

The **Execution Store** will keep track of statuses as the workflow runs:

| | `single_task` | `string_out` |
|---|---|---|
|1|`NotStarted`|`NotStarted`|
|2|`QueuedInCromwell`|`NotStarted`|
|3|`Starting`|`NotStarted`|
|4|`Running`|`NotStarted`|
|5|`Done`|`NotStarted`|
|6|`Done`|`Running`|
|7|`Done`|`Done`|

In step 1, the workflow has just started and the ExecutionStore is created in its initial
state. The Value Store doesn't track statuses and so begins empty: `{ }`.

In steps 2-4, the Execution Store tracks the `single_task` job as the engine is executing it.

As the Execution Store is updated to indicate task completion is step 5, the Value Store is also updated to
include the output value of the task:
```json
{
  "single_task.string_out": "hello"
}
```

By step 6, Cromwell can use the fact that the task is complete to decide that the output node is ready to be
evaluated. And the input to the output expression is available for lookup in the Value Store.

In step 7, all workflow nodes have run and the workflow is complete. The Value Store is updated once again to
additionally contain the output node value:
```json
{
  "single_task.string_out": "hello",
  "string_out": "hello"
}
``` 

Cromwell can use this information to trigger the "workflow complete" logic.

## Handling scatters

When Cromwell runs scattered tasks, the Execution Store cannot tell ahead of time how many
`JobKey`s it will need to represent all of the shards in the scatter. It can get around
this problem by putting a placeholder `JobKey` for the scatter node in the Execution Store. When
the scatter key is evaluated, it expands the Execution Store to include new `JobKey`s representing
every shard in the scatter.

As with the single task example, the Value Store starts empty, and is updated with the results of each
shard only as and when they are generated.

To see that in action, Consider this workflow:

```wdl
version 1.0

workflow scattered_task_workflow {
  scatter (x in range(2)) {
    call scattered_task
  }
  output {
    Int results_count = length(scattered_task.string_out)
  }
}

task scattered_task {
  command {
    echo hello
  }
  output {
    String string_out = "hello"
  }
}
```

#### Scatter Expansion

As the workflow starts, the execution store has three entries. An `x` represents the array-input for the scatter,
a`ScatterNode` represents the placeholder for expanding the scatter, and a `results_count` represents the workflow output. 

The start of workflow execution looks like this:

| | `x` | `ScatterNode` | `results_count` |
|---|---|---|---|
|1|`NotStarted`|`NotStarted`| `NotStarted` |
|2|`Running`|`NotStarted`| `NotStarted` | 
|3|`Done`|`NotStarted`| `NotStarted` | 

Once `x` is evaluated the value store gains an entry:
```json
{
  "x": "[0, 1]"
}
```

The scatter node now becomes runnable because its upstream dependency (`x`) is `Done` in the Execution Store.

The evaluation of `ScatterNode` updates the execution store in a number of ways:

* One call key for each index of `scattered_task` is added.
* The `scattered_task` gets an un-indexed key too. This key is used to mark when all of the shards of the call are complete.
* The gathered value `scattered_task.string_out` represents the "gathered" results of the task's output. It only runs 
once the un-indexed `scattered_task` key is Done and gathers output values into an array.
This gather key also acts as the upstream dependency of the `results_count` output expression.
* The `ScatterNode` is marked as `Done` so that it doesn't get triggered to run again.

Following the scatter-expansion evaluation of `ScatterNode`, the Execution Store looks like this:

| | `x` | `ScatterNode` | `scattered_task:0` | `scattered_task:1` | `scattered_task` | `scattered_task.string_out` | `results_count` |
|---|---|---|---|---|---|---|---|
|4|`Done`|`Done`|`NotStarted`|`NotStarted`|`NotStarted`|`NotStarted`|`NotStarted`

The Value Store is not changed at this time because no new values have been generated.

#### Parallel Shard Execution

The two scattered shards are now immediately runnable because they have no upsteam dependencies.
As the two jobs are run, the Execution Store map updates to track their statuses:

| | `x` | `ScatterNode` | `scattered_task:0` | `scattered_task:1` | `scattered_task` | `scattered_task.string_out` | `results_count` |
|---|---|---|---|---|---|---|---|
|4|`Done`|`Done`|`NotStarted`|`NotStarted`|`NotStarted`|`NotStarted`|`NotStarted`
|5|`Done`|`Done`|`QueuedInCromwell`|`NotStarted`|`NotStarted`|`NotStarted`|`NotStarted`
|6|`Done`|`Done`|`QueuedInCromwell`|`QueuedInCromwell`|`NotStarted`|`NotStarted`|`NotStarted`
|7|`Done`|`Done`|`Starting`|`QueuedInCromwell`|`NotStarted`|`NotStarted`|`NotStarted`
|8|`Done`|`Done`|`Starting`|`Starting`|`NotStarted`|`NotStarted`|`NotStarted`
|9|`Done`|`Done`|`Running`|`Starting`|`NotStarted`|`NotStarted`|`NotStarted`
|10|`Done`|`Done`|`Running`|`Running`|`NotStarted`|`NotStarted`|`NotStarted`
|11|`Done`|`Done`|`Running`|`Done`|`NotStarted`|`NotStarted`|`NotStarted`
|12|`Done`|`Done`|`Done`|`Done`|`NotStarted`|`NotStarted`|`NotStarted`

As the results for each shard come in, the value store is also updated to include them:

At step 11 (shard 1 has finished but shard 0 has not):
```json
{
  "x": "[0, 1]",
  "scattered_task.string_out:1": "hello"
}
```

At step 12:
```json
{
  "x": "[0, 1]",
  "scattered_task.string_out:0": "hello",
  "scattered_task.string_out:1": "hello"
}
```

#### Scatter Completion and Gathering

Once all of the sharded keys for `scattered_task` are complete, the un-indexed marker key for that call becomes
runnable. And once the marker is complete, the gather key for the output also becomes runnable.

The progression in the Execution Store goes like:

| | `x` | `ScatterNode` | `scattered_task:0` | `scattered_task:1` | `scattered_task` | `scattered_task.string_out` | `results_count` |
|---|---|---|---|---|---|---|---|
|12|`Done`|`Done`|`Done`|`Done`|`NotStarted`|`NotStarted`|`NotStarted`
|13|`Done`|`Done`|`Done`|`Done`|`Done`|`NotStarted`|`NotStarted`
|14|`Done`|`Done`|`Done`|`Done`|`Done`|`Done`|`NotStarted`

As the gather node completes in step 14, the value store is also updated to contain the unindexed, gathered result of 
the `scattered_task.string_out` output:

```json
{
  "x": "[0, 1]",
  "scattered_task.string_out:0": "hello",
  "scattered_task.string_out:1": "hello",
  "scattered_task.string_out": ["hello", "hello"]
}
```

When the `scattered_task.string_out` gather node completes, the upstream dependencies of the `results_count` output are
finally satisfied and it becomes runnable too. It runs to produce the workflow outputs:

| | `x` | `ScatterNode` | `scattered_task:0` | `scattered_task:1` | `scattered_task` | `scattered_task.string_out` | `results_count` |
|---|---|---|---|---|---|---|---|
|14|`Done`|`Done`|`Done`|`Done`|`Done`|`Done`|`NotStarted`
|15|`Done`|`Done`|`Done`|`Done`|`Done`|`Done`|`Running`
|16|`Done`|`Done`|`Done`|`Done`|`Done`|`Done`|`Done`

For step 16, completion of the output evaluation creates an entry in the Value Store which can be exposed as a workflow output as the 
workflow completes:

```json
{
  "x": "[0, 1]",
  "scattered_task.string_out:0": "hello",
  "scattered_task.string_out:1": "hello",
  "scattered_task.string_out": ["hello", "hello"],
  "results_count": 2
}
```