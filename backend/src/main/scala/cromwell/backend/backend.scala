package cromwell.backend

import com.typesafe.config.Config
import cromwell.core.WorkflowOptions.WorkflowOption
import cromwell.core.callcaching.MaybeCallCachingEligible
import cromwell.core.labels.Labels
import cromwell.core.{CallKey, WorkflowId, WorkflowOptions}
import cromwell.services.keyvalue.KeyValueServiceActor.KvResponse
import wdl4s.wdl._
import wdl4s.wdl.values.WdlValue

import scala.util.Try

/**
  * For uniquely identifying a job which has been or will be sent to the backend.
  */
case class BackendJobDescriptorKey(call: WdlTaskCall, index: Option[Int], attempt: Int) extends CallKey {
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
                                inputDeclarations: EvaluatedTaskInputs,
                                maybeCallCachingEligible: MaybeCallCachingEligible,
                                prefetchedKvStoreEntries: Map[String, KvResponse]) {
  val fullyQualifiedInputs = inputDeclarations map { case (declaration, value) => declaration.fullyQualifiedName -> value }
  val call = key.call
  override val toString = s"${key.mkTag(workflowDescriptor.id)}"
}

object BackendWorkflowDescriptor {
  def apply(id: WorkflowId,
            workflow: WdlWorkflow,
            knownValues: Map[FullyQualifiedName, WdlValue],
            workflowOptions: WorkflowOptions,
            customLabels: Labels) = {
    new BackendWorkflowDescriptor(id, workflow, knownValues, workflowOptions, customLabels, List.empty)
  }
}

/**
  * For passing to a BackendActor construction time
  */
case class BackendWorkflowDescriptor(id: WorkflowId,
                                     workflow: WdlWorkflow,
                                     knownValues: Map[FullyQualifiedName, WdlValue],
                                     workflowOptions: WorkflowOptions,
                                     customLabels: Labels,
                                     breadCrumbs: List[BackendJobBreadCrumb]) {
  
  val rootWorkflow = breadCrumbs.headOption.map(_.workflow).getOrElse(workflow)
  val rootWorkflowId = breadCrumbs.headOption.map(_.id).getOrElse(id)
  
  override def toString: String = s"[BackendWorkflowDescriptor id=${id.shortString} workflowName=${workflow.unqualifiedName}]"
  def getWorkflowOption(key: WorkflowOption) = workflowOptions.get(key).toOption
}

/**
  * For passing to a BackendActor construction time
  */
case class BackendConfigurationDescriptor(backendConfig: Config, globalConfig: Config) {

  lazy val backendRuntimeConfig = if (backendConfig.hasPath("default-runtime-attributes"))
    Option(backendConfig.getConfig("default-runtime-attributes")) else None
}

final case class AttemptedLookupResult(name: String, value: Try[WdlValue]) {
  def toPair = name -> value
}
