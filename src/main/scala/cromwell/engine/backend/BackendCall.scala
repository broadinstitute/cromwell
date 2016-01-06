package cromwell.engine.backend

import akka.event.LoggingAdapter
import cromwell.engine._
import Hashing._
import cromwell.binding._
import cromwell.binding.expression.WdlStandardLibraryFunctions
import cromwell.binding.values.WdlValue
import cromwell.engine.workflow.CallKey
import cromwell.logging.WorkflowLogger

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

/**
 * Represents a Call that is intended to be run on a specific Backend.
 *
 * A BackendCall is created by calling `Backend.bindCall(c: Call)`.  The
 * combination of a Call and a Backend give enough context to be able to
 * instantiate the Call's Task's command line, and provide an implementation
 * of the WDL standard library functions as well as a lookup functions to
 * resolve identifiers in expressions.
 *
 * The BackendCall also includes a Map[String, WdlValue] (CallInputs) which
 * represents the locally-qualified input values for this call.
 *
 * BackendCall aims to package everything that a Call needs to run into one
 * object so a client can do a parameterless `.execute()` to run the job.
 *
 * To understand why one needs a Backend and locally-qualified inputs in order
 * to instantiate a command, consider this example:
 *
 * task test {
 *   Array[String] something
 *   command {
 *     ./analysis --something-list=${write_lines(something)}
 *   }
 * }
 *
 * Without the Backend and the locallyQualifiedInputs, the expression
 * `write_lines(something)` could not be evaluated because:
 *
 * 1) We need to resolve "something" in locallyQualifiedInputs
 * 2) We need a place to write the results of `write_lines()`, and that location
 *    is backend-specific
 */

/**
 * Trait to encapsulate whether an execution is complete and if so provide a result.  Useful in conjunction
 * with the `poll` API to feed results of previous job status queries forward.
 */
trait ExecutionHandle {
  def isDone: Boolean
  def result: ExecutionResult
}

final case class CompletedExecutionHandle(override val result: ExecutionResult) extends ExecutionHandle {
  override val isDone = true
}

final case class SuccessfulExecutionHandle(outputs: CallOutputs, events: Seq[ExecutionEventEntry], returnCode: Int, hash: ExecutionHash, resultsClonedFrom: Option[BackendCall] = None) extends ExecutionHandle {
  override val isDone = true
  override val result = SuccessfulExecution(outputs, events, returnCode, hash, resultsClonedFrom)
}

final case class FailedExecutionHandle(throwable: Throwable, returnCode: Option[Int] = None) extends ExecutionHandle {
  override val isDone = true
  override val result = FailedExecution(throwable, returnCode)
}

trait BackendCall {

  /**
   * The Workflow and Call to invoke.  It is assumed that in the creation
   * of a BackendCall object that the 'call' would be within the workflow
   */
  def workflowDescriptor: WorkflowDescriptor
  def key: CallKey
  def call = key.scope

  /**
   * Backend which will be used to execute the Call
   */
  def backend: Backend

  /**
   * Inputs to the call.  For example, if a call's task specifies a command like this:
   *
   * command {
   *   File some_dir
   *   ls ${some_dir}
   * }
   *
   * Then locallyQualifiedInputs could be Map("some_dir" -> WdlFile("/some/path"))
   */
  def locallyQualifiedInputs: CallInputs

  /**
   * Implementation of the WDL Standard Library Functions, used to evaluate WdlExpressions
   */
  def engineFunctions: WdlStandardLibraryFunctions

  /**
   * Function used to get the value for identifiers in expressions.  For example, the
   * expression `read_lines(my_file_var)` would have to call lookupFunction()("my_file_var")
   * during expression evaluation
   */
  def lookupFunction: String => WdlValue = WdlExpression.standardLookupFunction(locallyQualifiedInputs, key.scope.task.declarations, engineFunctions)

  /**
   * Attempt to evaluate all the ${...} tags in a command and return a String representation
   * of the command.  This could fail for a variety of reasons related to expression evaluation
   * which is why it returns a Try[String]
   */
  def instantiateCommand: Try[String] = {
    val backendInputs = backend.adjustInputPaths(key, locallyQualifiedInputs, workflowDescriptor)
    call.instantiateCommandLine(backendInputs, engineFunctions)
  }

  /** Initiate execution, callers can invoke `poll` once this `Future` completes successfully. */
  def execute(implicit ec: ExecutionContext): Future[ExecutionHandle]

  /**
   * The default implementation of this method is not expected to be called and simply throws an `NotImplementedError`.
   * If the corresponding backend does not override `Backend#findResumableExecutions` to return resumable executions,
   * this method will not be called.  If the backend does override `Backend#findResumableExecutions`, the corresponding
   * `BackendCall` should override this method to actually do the resumption work.
   */
  def resume(jobKey: JobKey)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
    throw new NotImplementedError(s"resume() called on a non-resumable BackendCall: $this")
  }

  def useCachedCall(cachedBackendCall: BackendCall)(implicit ec: ExecutionContext): Future[ExecutionHandle] = ???

  /** Given the specified value for the Docker hash, return the overall hash for this `BackendCall`. */
  private def hashGivenDockerHash(dockerHash: Option[String]): ExecutionHash = {
    val runtime = call.task.runtimeAttributes

    val orderedInputs = locallyQualifiedInputs.toSeq.sortBy(_._1)
    val orderedOutputs = call.task.outputs.sortWith((l, r) => l.name > r.name)
    val orderedRuntime = Seq(
      ("docker", dockerHash getOrElse ""),
      ("defaultZones", runtime.defaultZones.sorted.mkString(",")),
      ("failOnStderr", runtime.failOnStderr.toString),
      ("continueOnReturnCode", runtime.continueOnReturnCode match {
        case ContinueOnReturnCodeFlag(bool) => bool.toString
        case ContinueOnReturnCodeSet(codes) => codes.toList.sorted.mkString(",")
      }),
      ("cpu", runtime.cpu.toString),
      ("preemptible", runtime.preemptible.toString),
      ("defaultDisks", runtime.defaultDisks.sortWith((l, r) => l.getName > r.getName).map(d => s"${d.getName} ${d.size} ${d.getType}").mkString(",")),
      ("memoryGB", runtime.memoryGB.toString)
    )

    val overallHash = Seq(
      backend.backendType.toString,
      call.task.commandTemplateString,
      orderedInputs map { case (k, v) => s"$k=${v.computeHash(workflowDescriptor.fileHasher).value}" } mkString "\n",
      orderedRuntime map { case (k, v) => s"$k=$v" } mkString "\n",
      orderedOutputs map { o => s"${o.wdlType.toWdlString} ${o.name} = ${o.expression.toWdlString}" } mkString "\n"
    ).mkString("\n---\n").md5Sum

    ExecutionHash(overallHash, dockerHash)
  }

  /**
    * Compute a hash that uniquely identifies this call
    */
  def hash(implicit ec: ExecutionContext): Future[ExecutionHash] = {
    // If a Docker image is defined in the task's runtime attributes, return a `Future[Option[String]]` of the Docker
    // hash string, otherwise return a `Future.successful` of `None`.
    def hashDockerImage(dockerImage: String): Future[Option[String]] =
      if (workflowDescriptor.lookupDockerHash)
        backend.dockerHashClient.getDockerHash(dockerImage) map { dh => Option(dh.hashString) }
      else
        Future.successful(Option(dockerImage))

    if (workflowDescriptor.configCallCaching)
      call.task.runtimeAttributes.docker map hashDockerImage getOrElse Future.successful(None) map hashGivenDockerHash
    else
      Future.successful(ExecutionHash("", None))
  }

  /**
   * Using the execution handle from the previous execution, resumption, or polling attempt, poll the execution
   * of this `BackendCall`.
   */
  def poll(previous: ExecutionHandle)(implicit ec: ExecutionContext): Future[ExecutionHandle]

  def workflowLoggerWithCall(source: String, akkaLogger: Option[LoggingAdapter] = None) = WorkflowLogger(
    source,
    workflowDescriptor,
    akkaLogger = akkaLogger,
    callTag = Option(key.tag)
  )
}
