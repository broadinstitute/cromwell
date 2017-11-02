In both Run and Server mode, you have the ability to ask Cromwell to abort a running workflow. This section explains what that entails.

When aborting workflow, either through the [abort endpoint](api/RESTAPI#abort-a-running-workflow) or by terminating the [Cromwell run process](Modes) (if [configured](Configuring#abort) to do so), Cromwell will take the following steps:

- Change the status of the workflow to `Aborting`
- Do not start any new job going forward
- Ask every running job to abort
- Wait for all running jobs to reach completion
- Finalize the workflow
- Change the status of the workflow to `Aborted`

The action of aborting a job is backend specific. Cromwell can only ask a backend to abort a job and wait for the backend to notify when it is in fact aborted.
Note that by the time the backend is asked to abort a job, it might already be successful, or failed, in which case Cromwell will report the job as such.

As an example, the Google backend when asked to abort a job will send an abort request to the Pipelines API. When Pipelines API indicates that the status of the job is aborted, Cromwell will mark it as such.
However the abort action is entirely dependent on the backend, in this particular case, Pipelines API does not guarantee the success of an abort request (see [Pipelines API documentation on abort](https://cloud.google.com/genomics/reference/rest/v1alpha2/operations/cancel)).

You'll also notice that the workflow is finalized even though being aborted.
Finalization is the last step in the execution of a workflow and a chance for each backend to do some work before the workflow is terminated.
Abort won't deny the chance for backends to finalize the workflow.

## Restart

If Cromwell is running in server mode and is restarted while a workflow was `Aborting`, Cromwell will attempt to reconnect to all the jobs that were being aborted in that workflow such that their status can be updated.
The ability to reconnect to an existing job is again backend dependent. The Google Backend and the HPC backends currently support it.
In this particular case, because there is no guarantee that all running jobs had a chance to abort before Cromwell was restarted, the backends will be asked to reconnect and abort the job.

During this reconnection process, Cromwell might ask backends to reconnect to jobs that were never started before the restart. In that case, the job will be mark as failed with an explanation message. This failure is benign and only an artifact of the fact that Cromwell was restarted while the workflow was aborting.

If the backend does not support reconnection to an existing job, jobs will be marked as failed with an explanation message as well. The backend status of the jobs will be "Unknown".
