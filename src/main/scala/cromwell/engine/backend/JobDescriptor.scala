package cromwell.engine.backend

import cromwell.engine.WorkflowDescriptor
import cromwell.engine.workflow.{BackendCallKey, CallKey, FinalCallKey}
import cromwell.webservice.WorkflowMetadataResponse
import wdl4s._
import wdl4s.values.WdlValue


/** Aspires to be an equivalent of `TaskDescriptor` in the pluggable backends world, describing a job in a way
  * that is complete enough for it to be executed on any backend and free of references to engine types.
  * Currently a ways to go in freedom from engine types.
  *
  * @tparam K CallKey subtype
  */
sealed trait JobDescriptor[K <: CallKey] {
  def workflowDescriptor: WorkflowDescriptor
  def key: K
  def locallyQualifiedInputs: CallInputs
}

final case class BackendCallJobDescriptor(workflowDescriptor: WorkflowDescriptor,
                                          key: BackendCallKey,
                                          locallyQualifiedInputs: CallInputs = Map.empty) extends JobDescriptor[BackendCallKey] {

  // PBE temporarily still required.  Once we have call-scoped Backend actors they will know themselves and the
  // backend won't need to be in the WorkflowDescriptor and this method won't need to exist.
  def backend = workflowDescriptor.backend

  def callRootPath = backend.callRootPath(this)

  def callRootPathWithBaseRoot(baseRoot: String) = backend.callRootPathWithBaseRoot(this, baseRoot)
}

final case class FinalCallJobDescriptor(workflowDescriptor: WorkflowDescriptor,
                                        key: FinalCallKey,
                                        workflowMetadataResponse: WorkflowMetadataResponse) extends JobDescriptor[FinalCallKey] {

  override val locallyQualifiedInputs = Map.empty[String, WdlValue]
}

