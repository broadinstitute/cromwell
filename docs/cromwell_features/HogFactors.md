# Hog Factors

## Introduction

Cromwell has only a finite amount of resources at its disposal. WDL and CWL workflows allow scattered tasks to be
run a huge number of times with very simple syntax.
This makes it easy for a very small number of workflows to hog all of the resources of Cromwell, forcing all other 
workflows (even a simple 'hello_world') to wait in line behind them. 
Sometimes that's exactly what you want, but often you would like Cromwell to remain responsive to new users' small 
workflows even while continuing to process production workflows from established users.

Cromwell 35 provides new ways of stopping one workflow or group of workflows from locking out everyone else by 
introducing hog factors, hog groups and hog limits. This page describes what they are and how they work.

## Concepts

### Hog Group 

The Hog Group is a way of grouping workflow from different submissions together to restrain their overall resource 
usage as a whole.

- Every top-level workflow is assigned to a hog group when Cromwell receives it. 
    + Exactly how this happens is [configurable](#configuration).
- Every sub-workflow or call started by a workflow is associated with the same hog group as its parent workflow.
- Multiple top-level workflows can be assigned to the same hog group to allow them to be grouped together.

Thus:

- Every workflow is assigned to a hog group when submitted.
- Every hog group may have many workflows assigned to it.

### Hog Factor

The hog factor is an integer greater than or equal to 1. It represents a trade-off between: 

- Fully utilizing all resources available to Cromwell to complete jobs, for as long as there are jobs to be processed.
- Reserving resources for requests from other hog groups - even if we have jobs waiting to run that could be using them.

Here are a few mental models which might be helpful to thinking about the hog factor:

- A hog factor of 2 means that "2 greedy users would be able to hog the entire resources of Cromwell" 
- A hog factor of 100 means "any 1 group is only ever allowed to use 1/100th of the resources of the total Cromwell server"

### Hog Limit

A Hog Limit is how much of a given resource a hog group is allowed use. Hog limits are not set directly; they are 
values that Cromwell calculates internally.
For example, a single hog group may be limited by Cromwell to a hog limit of 200 jobs per group. Therefore no matter how 
many workflows, sub-workflows, and jobs are queued in a greedy hog group, the whole group is limited to 200 concurrent 
running jobs.


## Configuration

Cromwell accepts the following configuration values for hog factors in the `hog-safety` stanza of `system` in the configuration 
file:
```conf
system {
  hog-safety {
    hog-factor = 1
    workflow-option = "hogGroup"
    token-log-interval-seconds = 0
  }
}
```

### Setting a hog-factor

The hog factor option sets the integer described in the [Hog Factor](#hog-factor) section above.
The default value is `1` (which is equivalent to not limiting by hog group). 

### Assignment of hog groups

Within the configuration file, you can specify the workflow option that will determine the hog group. The default is
`hogGroup`. So if a workflow arrives with the following workflow options file, Cromwell will assign `hogGroupA` as the 
workflow's hog group:

```json
{
  "hogGroup": "hogGroupA"
}
```

- Any workflow option value can be used so long as it is a simple `String` value:
    + You can come up with a new field and set it specifically for assigning hog groups.
    + You can choose a field that is already being used for other reasons
- If a workflow is submitted without a value for the designated field in its workflow options, the workflow ID is used 
as the hog group identifier. 

### Logging

Because the system is not a simple first-in-first-out, it can be valuable to see the status of all 
the existing queues inside Cromwell. 

To have this information logged on a regular basis, you can enable periodic queue logging. This will
also alert you at the same frequency when events are happening, such as hog groups being at their individual limits.

* You can enable logging for the Job Execution Token Dispenser, using
the `system.hog-safety.token-log-interval-seconds` configuration value.
* The default, `0`, means that no logging will occur.

## Effects

### Job Execution

#### Reserving Space

- Cromwell allows administrators to designate an overall maximum concurrent job limit per
 [backend](../backends/Backends.md#backend-job-limits). 
- Within that limit, a hog factor allows us to limit the maximum concurrent jobs started *per hog group*.
    + This is what allows new jobs to run immediately even if many jobs from bigger workflows are already queued up.

#### Round robin allocation

Rather than starting jobs on a strict first-come first served basis, Cromwell now assigns in a round-robin
fashion between hog groups and then on a first-come-first-served *within* a hog group.

In other words if the hog groups had the following entries queued up:
```
 A: jobA1, jobA2, jobA3, ..., jobA1000000
 B: jobB1, jobB2
 C: jobC1
 D: jobD1, jobD2
```

Then Cromwell would start the jobs in the following order, even though `jobA1000000` was added before `jobD1`:
```
jobA1, jobB1, jobC1, jobD1, jobA2, jobB2, jobD2, jobA3, ..., jobA1000000
```

#### Example: How job execution is affected by hog factors

##### An administrator sets up a Cromwell server

- A Cromwell administrator sets the overall maximum concurrent job limit to 100,000 PAPIv2 jobs.
- The administrator also sets the hog factor to be 25.
- Cromwell will therefore calculate a per-hog-group concurrent job limit of 4,000 PAPIv2 jobs.

##### Our first hog group hits its limit

- 100 workflows are running in hog group "A" and between them have generated 20,000 jobs for PAPIv2. 
    + Cromwell initially starts 4,000 jobs.
    + Cromwell then starts the remaining 16,000 new jobs as existing jobs from this group finish.
    + New workflows in this group will not be able to start jobs either
        + Their jobs are queued behind the existing jobs from this hog group.
    + Note that Cromwell is currently only using 1/25th of its overall limit because its hog factor is 25.

##### Another hog group appears

- Now hog group B submits 1,000 workflows and between them they generate 200,000 jobs for PAPIv2.
- Even though 16,000 of group A's jobs are still queued, Cromwell starts 4,000 of group B's jobs immediately,
- The remaining 196,000 of group B's jobs are queued up waiting for group B's existing jobs to complete.

##### Where do we stand?

- Cromwell knows about 220,000 jobs that could be started
- Cromwell has an overall limit of 100,000
- Cromwell is running 8,000 jobs in two hog groups.

*In other words, not so great - perhaps we should have set the hog factor lower...?*

#####But wait, more workflows appear...

- Now another 23 hog groups ("C" through "Y") submit workflows of a similar scale to hog group A.
- One by one, the workflows of each hog group fill up their share of the overall concurrent job limit.
- So Cromwell is now running 100,000 jobs and each hog group has been allocated 4,000 of those.

##### What about poor hog group "Z"?

- A final group submits workflows under hog group "Z".
- Alas, even though hog group "Z" is not running anything yet, we cannot start their workflows because we're 
at the global maximum of 100,000.


*In other words, perhaps we should have set the hog factor higher...?*

##### So what now?

- As jobs in other hog group complete, we will begin to see hog group "Z" jobs started alongside new jobs from the 
other hog groups.
- Going forward Cromwell will start jobs from all groups at the same rate, even though hog group Z's jobs arrived
later than those from hog group A. Thus, over time, each group will approach approximately 1/26th of the total pool.
  
## FAQs

#### Can I opt out of using hog groups?

Yes, to various degrees:

- No matter what, your workflows will be assigned to a hog group. 
- To opt out of reserving Cromwell's resources for new hog groups, leave the hog factor set to 1.
- To opt out of round-robin allocation between workflows, and preserve a strict first-in-first-out allocation of jobs,
assign all workflows to the same hog-group in their workflow options.
    + To set this as the default, you can add a value to the default workflow options. 
    + For an example see the `workflow-options` / `default` stanza of [cromwell.examples.conf][cromwell-examples-conf].

[cromwell-examples-conf]: https://www.github.com/broadinstitute/cromwell/tree/develop/cromwell.examples.conf
