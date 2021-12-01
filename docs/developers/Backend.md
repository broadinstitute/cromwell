## Cromwell Backend Development

### What’s a Backend?

Cromwell is a workflow execution service, which means it takes a representation of work to do (tasks) as well as the
dependencies and code flow between them. When it comes to actually executing the specific task, Cromwell delegates that
work to a “backend”. A backend is therefore responsible for the actual execution of a single task on some underlying
computing platform, such as Google Cloud Platform, Azure, AWS, TES, SLURM, etc. Backends are implemented as a software
layer between the Cromwell engine and that underlying platform.

The underlying platform services requests for a job to run while the backend shim provides the interface layer between
Cromwell and the platform. In general, the more sophisticated the platform the thinner the shim needs to be and vice
versa.

A job from Cromwell’s perspective is a collection of the following information:

* A unix command line to run
* A mapping of where its input files currently live to where the command line expects them to be
* A mapping of where the command line will write its outputs to where they should eventually wind up
* An optional Docker image which if supplied will be the environment in which the command will be run
* A collection of arbitrary key/value pairs which can be used by the execution platform to tune the request, e.g. amount
  of memory or the number of CPUs.

The Docker image is not required for a backend but is highly recommended. The bioinformatics workflow field is rapidly
moving towards a Docker model so one will find better support for their backend if it is built around using Docker
containers.

The input/output file mappings will depend on the needs of the platform. In the simplest example these values would be
identical, living on a shared filesystem. However an example of where this would be useful would be a model where an
input file lives in a cloud bucket and is directly copied onto the machine running the command line followed by copying
the program’s outputs back to a cloud bucket.

The underlying platform can be extremely simple or have arbitrary levels of complexity as long as one can map from the
above concepts to running a command. That mapping could happen completely in the platform, the Cromwell backend, or a
mixture of the two. The implementer of a backend will need to strike the balance which works best for both their needs
and the particulars of the platform itself.

Throughout the rest of this document the following three terms are used, the difference between them is subtle but
important:

* Task: The abstract definition of a thing to run. Think of this like a function in a programming language.
* Call: An instantiated request in Cromwell of a thing to run. To further the function analogy this would be an
  invocation of that function.
* Job: The physical manifestation of a call by the backend. Examples of this would be an SGE job, a unix process, or a
  Google Life Sciences operation ID.

### Backend Lifecycle

Within a running workflow, the backend is used in the following manner:

* Initialization: Initialization routine called the first time a workflow uses a backend
* Execute: The workflow requests that the backend run a job
* Recover: Attempt to reconnect to a previously started job. An example of this would be if Cromwell was restarted and
  wanted to reattach to currently running jobs.
* Abort: Request that a running job be halted
* Finalization: When a workflow is complete, allows for any workflow level cleanup to take place.

The initialization and finalization steps are optional and are provided for cases where a backend needs that behavior;
not all backends will need these.

The implementations of both the recover and abort steps are up to the backend developer. For instance the recover
function could be implemented to actually perform an execution and/or the abort function could be written to do nothing
at all. Implementing these in a more robust manner is recommended for most cases but neither are universally
appropriate.

### How do I create a Backend?

This section is assuming that you both know the underlying execution platform you wish to use and that you know how to
programmatically submit work to it.

The Cromwell engine uses three different types of Akka actors in the backend while processing a workflow:

* `BackendWorkflowInitializationActor`: Handles the initialization phase
* `BackendJobExecutionActor`: Handles requests to execute a new job or recover an existing job
* `BackendWorkflowFinalizationActor`: Handles the finalization phase

These three actors are all represented by a trait which a backend developer would need to implement. Both the
initialization and finalization actors are optional.

The following explanations of the traits mention multiple types defined in the Cromwell codebase. In particular
`BackendJobDescriptor` wraps all the information necessary for a backend to instantiate a job and `JobKey` provides the
information to uniquely identify a job. It is recommended that one look in the Cromwell codebase for more information on
these and other types.

#### BackendWorkflowInitializationActor

If a backend developer wishes to take advantage of the initialization phase of the backend lifecycle they must implement
this trait. There are three functions which must be implemented:

1. `abortInitialization: Future[WorkflowAbortResponse]` specifies what to do, if anything, when a workflow is requested
   to abort while a backend initialization is in progress.
2. `validate: Future[Unit]` is provided so that the backend can ensure that all of the calls it will handle conform to
   the rules of that backend. For instance, if a backend requires particular runtime attributes to exist.
3. `beforeAll: Future[Unit]` is the actual initialization functionality. If a backend requires any work to be done prior
   to handling a call, that code must be called from here.

#### BackendJobExecutionActor

There are three functions to override in this trait

* `execute(jobDescriptor: BackendJobDescriptor): Future[BackendJobExecutionResponse]` requests that the backend run a
  job. This will require detailed knowledge of how the execution platform works and the implementation of this function
  will likely be unique to every platform. For instance an SGE backend using a shared filesystem might not need to do
  anything with filename mappings and simply submit the job’s command to the SGE server. Another example would be the
  [Google Life Sciences API](https://cloud.google.com/life-sciences/docs/reference/rest) which will spin up a VM with an
  instance type determined from specified runtime attributes (# of CPU, amount of RAM), pull the task’s Docker image,
  download the input files from Google Cloud Storage onto the VM, run the command in the Docker, copy the output files
  up to Google Cloud Storage and shut down the VM.
* `recover(jobDescriptor: BackendJobDescriptor): Future[BackendJobExecutionResponse]` requests that the backend
  reconnect to a previously launched job. It is recommended that a backend developer implement this as requested but
  that behavior isn’t appropriate for all situations, in which case it could for instance just call execute. In the case
  of an SGE backend this might be implemented to determine the SGE ID for the job and start from there.
* `abortJob(jobKey: BackendJobDescriptorKey): Future[JobAbortResponse]` requests that a running job be aborted. It is
  recommended that this be properly implemented but that’s not always appropriate. This function could easily be
  implemented as a no-op.

#### BackendWorkflowFinalizationActor

There are two functions to override for this trait if the backend developer chooses to use the finalization
functionality

* `abortFinalization: Future[WorkflowAbortResponse]` is called if a workflow is requested to abort while finalization is
  happening.
* `afterAll: Future[Unit]` is the hook to perform any desired functionality, and is called when the workflow is
  completed.
