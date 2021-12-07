## Cromwell Backend Development

### What’s a Backend?

Cromwell is a workflow execution service, which means it takes a representation of work to do (tasks) as well as the
dependencies and code flow between them. When it comes to actually executing the specific task, Cromwell delegates that
work to a “backend”. A backend is therefore responsible for the actual execution of a single task on some underlying
computing platform, such as Google Cloud Platform, AWS, TES, Local, etc. Backends are implemented as a software layer
between the Cromwell engine and that underlying platform.

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

Throughout the rest of this document the following terms are used, the difference between them may be subtle but
important:

* Workflow: a container of one or more calls, possibly expressing execution dependencies on each other. May contain
  scatters, conditional logic, or subworkflow invocations.
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

The Cromwell engine uses backend-specific types extending Akka `Actor` while processing a workflow:

* [`BackendLifecycleActorFactory`](https://github.com/broadinstitute/cromwell/blob/9bf1622ca8988365477b77b9f26ce388b54fc58c/backend/src/main/scala/cromwell/backend/BackendLifecycleActorFactory.scala#L17)
  The entry point into the backend implementation specified in Cromwell configuration. e.g.
  a [sample configuration for the Local backend](https://github.com/broadinstitute/cromwell/blob/2b19f00976ee258142185917083460d724f7fe3d/cromwell.example.backends/cromwell.examples.conf#L370)
* [`BackendWorkflowInitializationActor`](https://github.com/broadinstitute/cromwell/blob/93392acf2881921dcf22ef4dbda12af42339b3ab/backend/src/main/scala/cromwell/backend/BackendWorkflowInitializationActor.scala#L27)
  Handles the initialization phase
* Usually a derivative
  of [`StandardAsyncExecutionActor`](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L75)
  to run an individual job within the workflow
* [`BackendWorkflowFinalizationActor`](https://github.com/broadinstitute/cromwell/blob/a40de672c565c4bbd40f57ff96d4ee520dc2b4fc/backend/src/main/scala/cromwell/backend/BackendWorkflowFinalizationActor.scala#L10)
  Handles the finalization phase

These three actors are all represented by a trait which a backend developer would need to implement. Both the
initialization and finalization actors are optional.

The following explanations of the traits mention multiple types defined in the Cromwell codebase. In particular
`BackendJobDescriptor` wraps all the information necessary for a backend to instantiate a job and `JobKey` provides the
information to uniquely identify a job. It is recommended that one look in the Cromwell codebase for more information on
these and other types.

#### BackendWorkflowInitializationActor

If a backend developer wishes to take advantage of the initialization phase of the backend lifecycle they must implement
this trait. There are three functions which must be implemented:

* [`abortInitialization: Unit`](https://github.com/broadinstitute/cromwell/blob/93392acf2881921dcf22ef4dbda12af42339b3ab/backend/src/main/scala/cromwell/backend/BackendWorkflowInitializationActor.scala#L168)
   specifies what to do, if anything, when a workflow is requested to abort while a backend initialization is in
   progress.
* [`validate: Future[Unit]`](https://github.com/broadinstitute/cromwell/blob/93392acf2881921dcf22ef4dbda12af42339b3ab/backend/src/main/scala/cromwell/backend/BackendWorkflowInitializationActor.scala#L178)
   is provided so that the backend can ensure that all of the calls it will handle conform to the rules of that backend.
   For instance, if a backend requires particular runtime attributes to exist.
* [`beforeAll: Future[Option[BackendInitializationData]]`](https://github.com/broadinstitute/cromwell/blob/93392acf2881921dcf22ef4dbda12af42339b3ab/backend/src/main/scala/cromwell/backend/BackendWorkflowInitializationActor.scala#L173)
   is the actual initialization functionality. If a backend requires any work to be done prior to handling a call, that
   code must be called from here.

#### StandardAsyncExecutionActor

Nearly all production backend implementations in Cromwell extend
the [StandardAsyncExecutionActor](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L75)
trait. The minimum overrides for implementations of this trait are enumerated below, but overrides of other methods will
also be required. There are several backend implementations in the Cromwell codebase that can serve as references, for
example:

* [Google Life Sciences / Pipelines API](https://github.com/broadinstitute/cromwell/blob/0aff35336b4e2ba19b18530a68e622df1462d9b7/supportedBackends/google/pipelines/common/src/main/scala/cromwell/backend/google/pipelines/common/PipelinesApiAsyncBackendJobExecutionActor.scala#L95)
  common layer
    * [Life Sciences API (beta)](https://github.com/broadinstitute/cromwell/blob/a49e1fc65703ccfda2840d1d9266fad2bdbb7339/supportedBackends/google/pipelines/v2beta/src/main/scala/cromwell/backend/google/pipelines/v2beta/PipelinesApiAsyncBackendJobExecutionActor.scala#L27)
    * [Pipelines API (alpha)](https://github.com/broadinstitute/cromwell/blob/a49e1fc65703ccfda2840d1d9266fad2bdbb7339/supportedBackends/google/pipelines/v2alpha1/src/main/scala/cromwell/backend/google/pipelines/v2alpha1/PipelinesApiAsyncBackendJobExecutionActor.scala#L27)
* [GA4GH Task Execution Service (TES)](https://github.com/broadinstitute/cromwell/blob/6bf7af3c12a411db26786ac34646238fc053ec97/supportedBackends/tes/src/main/scala/cromwell/backend/impl/tes/TesAsyncBackendJobExecutionActor.scala#L55)
* [AWS Batch](https://github.com/broadinstitute/cromwell/blob/470d482e8ba2a9e2bc544896a4e6ceea57d55bb2/supportedBackends/aws/src/main/scala/cromwell/backend/impl/aws/AwsBatchAsyncBackendJobExecutionActor.scala#L83)

Overrides required for compilation and basic execution of `StandardAsyncExecutionActor` implementations:

* [`type StandardAsyncRunInfo`](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L89)
  encapsulates the type of the run info when a job is started.
* [`type StandardAsyncRunState`](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L92)
  encapsulates the type of the run status returned during each poll.
* [`def statusEquivalentTo(thiz: StandardAsyncRunState)(that: StandardAsyncRunState): Boolean`](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L97)
  should return true if the status contained in `thiz` is equivalent to `that`, delta any other data that might be
  carried around in the state type
* [`def standardParams: StandardAsyncExecutionActorParams`](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L103)
  a standard set of parameters passed to the backend.
* [`def isTerminal(runStatus: StandardAsyncRunState): Boolean`](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L829)
  Returns true when a job is complete, either successfully or unsuccessfully.
* At least one of:
  * [`def execute(): ExecutionHandle`](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L738) executes the job specified in the params
  * [`def executeAsync(): Future[ExecutionHandle]`](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L748) asynchronously executes the job specified in the params
* At least one of:
  * [`def pollStatus(handle: StandardAsyncPendingExecutionHandle): StandardAsyncRunState`](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L798) returns the run status for the job
  * [`def pollStatusAsync(handle: StandardAsyncPendingExecutionHandle): Future[StandardAsyncRunState]`](https://github.com/broadinstitute/cromwell/blob/9181235d364712b78dbea1f35042c3c6e431af87/backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala#L808) asynchronously returns the run status for the job

#### BackendWorkflowFinalizationActor

There is only one function to override for this trait if the backend developer chooses to use the finalization
functionality:

* [`afterAll: Future[Unit]`](https://github.com/broadinstitute/cromwell/blob/a40de672c565c4bbd40f57ff96d4ee520dc2b4fc/backend/src/main/scala/cromwell/backend/BackendWorkflowFinalizationActor.scala#L37)
  is the hook to perform any desired functionality, and is called when the workflow is completed.
