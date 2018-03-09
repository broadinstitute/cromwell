package cromwell.backend

import _root_.wdl.draft2.model._
import com.typesafe.config.Config
import cromwell.core.WorkflowOptions.WorkflowOption
import cromwell.core.callcaching.MaybeCallCachingEligible
import cromwell.core.labels.Labels
import cromwell.core.{CallKey, WorkflowId, WorkflowOptions}
import cromwell.services.keyvalue.KeyValueServiceActor.KvResponse
import wom.callable.ExecutableCallable
import wom.graph.GraphNodePort.OutputPort
import wom.graph.CommandCallNode
import wom.values.{WomEvaluatedCallInputs, WomValue}

import scala.util.Try

/**
  * For uniquely identifying a job which has been or will be sent to the backend.
  */
case class BackendJobDescriptorKey(call: CommandCallNode, index: Option[Int], attempt: Int) extends CallKey {
  def node = call
  private val indexString = index map { _.toString } getOrElse "NA"
  lazy val tag = s"${call.fullyQualifiedName}:$indexString:$attempt"
  def mkTag(workflowId: WorkflowId) = s"$workflowId:$this"
}

/**
  * For passing to a BackendWorkflowActor for job execution or recovery
  */
case class BackendJobDescriptor(workflowDescriptor: BackendWorkflowDescriptor,
                                key: BackendJobDescriptorKey,
                                runtimeAttributes: Map[LocallyQualifiedName, WomValue],
                                evaluatedTaskInputs: WomEvaluatedCallInputs,
                                maybeCallCachingEligible: MaybeCallCachingEligible,
                                prefetchedKvStoreEntries: Map[String, KvResponse]) {
  val fullyQualifiedInputs = evaluatedTaskInputs map { case (declaration, value) =>
    key.call.identifier.combine(declaration.name).fullyQualifiedName.value -> value
  }
  val localInputs = evaluatedTaskInputs map { case (declaration, value) => declaration.name -> value }
  val taskCall = key.call
  override lazy val toString = key.mkTag(workflowDescriptor.id)
}

object BackendWorkflowDescriptor {
  def apply(id: WorkflowId,
            callable: ExecutableCallable,
            knownValues: Map[OutputPort, WomValue],
            workflowOptions: WorkflowOptions,
            customLabels: Labels) = {
    new BackendWorkflowDescriptor(id, callable, knownValues, workflowOptions, customLabels, List.empty)
  }
}

/**
  * For passing to a BackendActor construction time
  */
case class BackendWorkflowDescriptor(id: WorkflowId,
                                     callable: ExecutableCallable,
                                     knownValues: Map[OutputPort, WomValue],
                                     workflowOptions: WorkflowOptions,
                                     customLabels: Labels,
                                     breadCrumbs: List[BackendJobBreadCrumb]) {

  val rootWorkflow = breadCrumbs.headOption.map(_.callable).getOrElse(callable)
  val rootWorkflowId = breadCrumbs.headOption.map(_.id).getOrElse(id)

  override def toString: String = s"[BackendWorkflowDescriptor id=${id.shortString} workflowName=${callable.name}]"
  def getWorkflowOption(key: WorkflowOption) = workflowOptions.get(key).toOption
}

/**
  * For passing to a BackendActor construction time
  */
case class BackendConfigurationDescriptor(backendConfig: Config, globalConfig: Config) {

  lazy val backendRuntimeConfig =
    if (backendConfig.hasPath("default-runtime-attributes"))
      Option(backendConfig.getConfig("default-runtime-attributes"))
    else
      None
}

final case class AttemptedLookupResult(name: String, value: Try[WomValue]) {
  def toPair = name -> value
}
