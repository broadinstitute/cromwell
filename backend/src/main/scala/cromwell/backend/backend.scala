package cromwell.backend

import com.typesafe.config.Config
import cromwell.core.WorkflowOptions.WorkflowOption
import cromwell.core.{CallKey, WorkflowId, WorkflowOptions}
import wdl4s._
import wdl4s.values.WdlValue

import scala.util.Try

/**
  * For uniquely identifying a job which has been or will be sent to the backend.
  */
case class BackendJobDescriptorKey(call: TaskCall, index: Option[Int], attempt: Int) extends CallKey {
  def scope = call
  private val indexString = index map { _.toString } getOrElse "NA"
  val tag = s"${call.fullyQualifiedName}:$indexString:$attempt"
  def mkTag(workflowId: WorkflowId) = s"$workflowId:$this"
}

/**
  * For passing to a BackendWorkflowActor for job execution or recovery
  */
case class BackendJobDescriptor(workflowDescriptor: BackendWorkflowDescriptor,
                                key: BackendJobDescriptorKey,
                                runtimeAttributes: Map[LocallyQualifiedName, WdlValue],
                                inputDeclarations: EvaluatedTaskInputs) {
  val fullyQualifiedInputs = inputDeclarations map { case (declaration, value) => declaration.fullyQualifiedName -> value }
  val call = key.call
  override val toString = s"${key.mkTag(workflowDescriptor.id)}"
}

object BackendWorkflowDescriptor {
  def apply(id: WorkflowId,
            workflow: Workflow,
            knownValues: Map[FullyQualifiedName, WdlValue],
            workflowOptions: WorkflowOptions) = {
    new BackendWorkflowDescriptor(id, workflow, knownValues, workflowOptions, List.empty)
  }
}

/**
  * For passing to a BackendActor construction time
  */
case class BackendWorkflowDescriptor(id: WorkflowId,
                                     workflow: Workflow,
                                     knownValues: Map[FullyQualifiedName, WdlValue],
                                     workflowOptions: WorkflowOptions,
                                     breadCrumbs: List[BackendJobBreadCrumb]) {
  
  val rootWorkflow = breadCrumbs.headOption.map(_.workflow).getOrElse(workflow)
  val rootWorkflowId = breadCrumbs.headOption.map(_.id).getOrElse(id)
  
  override def toString: String = s"[BackendWorkflowDescriptor id=${id.shortString} workflowName=${workflow.unqualifiedName}]"
  def getWorkflowOption(key: WorkflowOption) = workflowOptions.get(key).toOption
}

/**
  * For passing to a BackendActor construction time
  */
case class BackendConfigurationDescriptor(backendConfig: Config, globalConfig: Config)

final case class AttemptedLookupResult(name: String, value: Try[WdlValue]) {
  def toPair = name -> value
}

case class PreemptedException(msg: String) extends Exception(msg)
