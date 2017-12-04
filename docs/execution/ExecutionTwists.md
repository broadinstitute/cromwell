In order to run a workflow, Cromwell uses the backends available to it to create jobs and monitor them until they are complete. When no new jobs can start and all jobs have finished (Success, Failure, Aborted), the workflow terminates.

In an ideal situation, all jobs succeed and the workflow succeeds.
Unfortunately things don't always go as planned. Here is what to expect from Cromwell when things go off the rails.

## Failure Modes

Cromwell supports two failure modes, which specify how Cromwell behaves when a job fails during the execution of a workflow.

* `NoNewCalls` **(default)**  
	* Cromwell does not start any new call as soon as a job fails. Cromwell will still monitor the rest of the jobs until they complete (successfully or not).  
* `ContinueWhilePossible`  
	* Cromwell attempts to run as many jobs as possible until no more can be started. When all running jobs are complete, the workflow fails.

The failure mode can be set in the [Configuration](../Configuring/) or [Workflow Options](../wf_options/Overview#workflow-failure).

_For example:_

![](ABdependency.png)

This simple diagram represents 4 jobs:

* Job A and Job B are independent; 
* Job A1 depends on A; 
* Job B1 depends on B.

Let's look at the case where A and B are both running, and for some reason B fails (shaded red).

**`NoNewCalls`**  
If the failure mode is `NoNewCalls` Cromwell waits for A to complete (shaded green). Regardless of whether A is successful or not, Cromwell then fails the workflow, without starting A1 or B1 (shaded grey).

![](NNC_B_fail.png)

**`ContinueWhilePossible`**  
If the failure mode is `ContinueWhilePossible` and A succeeds (green), then Cromwell starts A1 (green) and waits for it to complete. At this point all jobs that can run have completed. B1 (grey) cannot run since B failed (red), therefore Cromwell fails the workflow without starting B1.

![](CWP_B_fail.png)

### Retryable failures

Retryable failures are **not** failures that can trigger a workflow failure. An example of a retryable failure is when a [preemptible VM](../RuntimeAttributes/#preemptible) is preempted. 

In the example above, if B's failure is retryable then B will be retried (shaded yellow and green). The workflow will keep running normally, regardless of which failure mode is enabled.

![](CWP_B_retryable_fail_then_success.png)

Using the previous example, let's imagine **B** failed from a **non-retryable** failure (shaded red). After B failed, **A** fails from a **retryable** failure. Now Cromwell's behavior will depend on the failure mode.

**`NoNewCalls`**  
If the failure mode is `NoNewCalls`, then A **will not be** retried (yellow). A1 (grey) will not start, because A did not complete.

![](NCC_B_fail_A_retryable.png)

**`ContinueWhilePossible`**  
If the failure mode is **`ContinueWhilePossible`**, then A **will be** retried (yellow and green). If A is successful then A1 will start (green).

![](CWP_B_fail_A_retryable.png)

## Abort

In both [Run](../Modes/#run) and [Server](../Modes/#server) mode, you can abort a running workflow. This section explains what that entails.

When aborting a workflow, either through the [abort endpoint](../api/RESTAPI#abort-a-running-workflow) or by terminating the [Cromwell run process](../Modes) (if [configured](../Configuring#abort) to do so), Cromwell does the following:

1. Changes the status of the workflow to `Aborting`,
2. Does not start any new jobs,
3. Asks every running job to abort,
4. Waits for all running jobs to complete,
5. Finalizes the workflow,
6. Changes the status of the workflow to `Aborted`.

The action of aborting a job is backend specific. Cromwell can only ask a backend to abort a job and wait for the backend to notify it when it is aborted.  

For example, if you are running Cromwell on the Google backend and abort a job, Google will send an abort request to the Google Pipelines API. When Pipelines API indicates that the status of the job is aborted, Cromwell will mark it as such.
Remember that abort action is entirely dependent on the backend. In this particular case, Pipelines API does not guarantee the success of an abort request (see [Pipelines API documentation on abort](https://cloud.google.com/genomics/reference/rest/v1alpha2/operations/cancel)).

You'll also notice that the workflow is finalized even though being aborted.
Finalization is the last step in the execution of a workflow and a chance for each backend to do some work before the workflow is terminated.
Backends won't be denied the chance to finalize the workflow even if it's being aborted.

_Note that by the time the backend is asked to abort a job, the job may have succeeded or failed already. In this case Cromwell will report the job's status (successful or failed)._

_If a job fails with a retryable failure (e.g is preempted), it will **not** be attempted again when the workflow is aborting._

## Restart

When Cromwell is run in server mode, it is likely that the server will have to be stopped and started again at some point (for example to upgrade Cromwell to a new version).
In such case, if the `workflow-restart` configuration option is set to `true` (default), Cromwell will restart all jobs that were in progress. Specifically Cromwell will attempt to reconnect to all running jobs of a workflow.

The ability to reconnect to an existing job is, again, backend dependent. The Google Backend and the HPC backends currently support it.

If the workflow was in state `Aborting`, Cromwell will ask all running jobs to abort again. No new jobs will be started.

Once all jobs have been reconnected to, the workflow will keep running normally.

During the reconnection process Cromwell might ask backends to reconnect to jobs that were never started before the restart. In that case, the job will be mark as failed with an explanation message. This failure is benign and only an artifact of the fact that Cromwell was restarted.
If the backend does not support reconnection to an existing job, jobs will be marked as failed with an explanation message as well. The backend status of the jobs will be "Unknown".
