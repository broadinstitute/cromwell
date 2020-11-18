# Workflow Execution: The Value Store

## Purpose

The Value Store is a data structure owned and maintained inside each 
`WorkflowExecutionActor` (see [Major Actors](majorActors.md)).

Its purpose is to hold the set of _values_ produced by the workflow so far. If
the WOM graph holds a static representation of the workflow which doesn't change
as the workflow is run, the Value Store records the values assigned to every 
task output and value definition evaluated so far during workflow execution.

**Note:** The Value Store does **not** hold the execution status of the various 
nodes in the WOM graph. Nor does it determine when downstream nodes are ready
to run. That is the domain of the [Execution Store](executionStore.md). 

## Data Structure

The Value Store data structure is a mapping from output ports on WOM `GraphNode`s 
(and shard index if necessary) to the appropriate `WomValue`. 

## Examples

Some worked through examples of how the Execution Store and Value Store change as workflows progress
are given on the [Execution and Value Store Examples](executionAndValueStoreExamples.md) page.