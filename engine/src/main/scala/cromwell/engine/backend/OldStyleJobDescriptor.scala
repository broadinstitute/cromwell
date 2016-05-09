package cromwell.engine.backend

import cromwell.backend.JobKey
import cromwell.CallEngineFunctions
import cromwell.engine.AbortRegistrationFunction
import cromwell.engine.backend.runtimeattributes.CromwellRuntimeAttributes
import cromwell.engine.workflow.{BackendCallKey, FinalCallKey}
import cromwell.webservice.WorkflowMetadataResponse
import wdl4s._
import wdl4s.values.WdlValue

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try


/** Aspires to be an equivalent of `TaskDescriptor` in the pluggable backends world, describing a job in a way
  * that is complete enough for it to be executed on any backend and free of references to engine types.
  * Currently a ways to go in freedom from engine types.
  *
  * @tparam K CallKey subtype
  */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
sealed trait OldStyleJobDescriptor[K <: JobKey] {
  def workflowDescriptor: OldStyleWorkflowDescriptor
  def key: K
  def locallyQualifiedInputs: CallInputs
}
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class OldStyleBackendCallJobDescriptor(workflowDescriptor: OldStyleWorkflowDescriptor,
                                            key: BackendCallKey,
                                            locallyQualifiedInputs: CallInputs = Map.empty,
                                            abortRegistrationFunction: Option[AbortRegistrationFunction] = None) extends OldStyleJobDescriptor[BackendCallKey] {

  lazy val call = key.scope

  // PBE temporarily still required.  Once we have call-scoped Backend actors they will know themselves and the
  // backend won't need to be in the WorkflowDescriptor and this method won't need to exist.
  lazy val backend = workflowDescriptor.backend

  lazy val callRootPath = backend.callRootPath(this)

  def callRootPathWithBaseRoot(baseRoot: String) = backend.callRootPathWithBaseRoot(this, baseRoot)

  lazy val callEngineFunctions: CallEngineFunctions = backend.callEngineFunctions(this)

  lazy val callRuntimeAttributes: CromwellRuntimeAttributes = backend.runtimeAttributes(this)

  def useCachedCall(cachedJobDescriptor: OldStyleBackendCallJobDescriptor)(implicit ec: ExecutionContext): Future[OldStyleExecutionHandle] = backend.useCachedCall(cachedJobDescriptor, this)

  /**
    * Attempt to evaluate all the ${...} tags in a command and return a String representation
    * of the command.  This could fail for a variety of reasons related to expression evaluation
    * which is why it returns a Try[String]
    */
  lazy val instantiateCommand: Try[String] = backend.instantiateCommand(this)

  def lookupFunction(evaluatedValues: Map[String, WdlValue]): String => WdlValue = {
    val currentlyKnownValues = locallyQualifiedInputs ++ evaluatedValues
    WdlExpression.standardLookupFunction(currentlyKnownValues, call.task.declarations, callEngineFunctions)
  }

  def poll(previous: OldStyleExecutionHandle)(implicit ec: ExecutionContext): Future[OldStyleExecutionHandle] = backend.poll(this, previous)

  /**
    * Compute a hash that uniquely identifies this job w.r.t. its backend.
    */
  def hash(implicit ec: ExecutionContext): Future[ExecutionHash] = backend.hash(this)

  def resume(executionInfos: Map[String, Option[String]])(implicit ec: ExecutionContext): Future[OldStyleExecutionHandle] = backend.resume(this, executionInfos)

  def execute(implicit ec: ExecutionContext): Future[OldStyleExecutionHandle] = backend.execute(this)
}
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class FinalCallJobDescriptor(workflowDescriptor: OldStyleWorkflowDescriptor,
                                  key: FinalCallKey,
                                  workflowMetadataResponse: WorkflowMetadataResponse) extends OldStyleJobDescriptor[FinalCallKey] {

  override val locallyQualifiedInputs = Map.empty[String, WdlValue]

  lazy val call = key.scope
}
