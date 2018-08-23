package cromwell.backend

import _root_.wdl.draft2.model._
import akka.actor.ActorSystem
import com.typesafe.config.Config
import common.validation.Validation._
import cromwell.core.WorkflowOptions.WorkflowOption
import cromwell.core.callcaching.MaybeCallCachingEligible
import cromwell.core.filesystem.CromwellFileSystems
import cromwell.core.labels.Labels
import cromwell.core.path.{DefaultPathBuilderFactory, PathBuilderFactory}
import cromwell.core.{CallKey, WorkflowId, WorkflowOptions}
import cromwell.services.keyvalue.KeyValueServiceActor.KvResponse
import wom.callable.{ExecutableCallable, MetaValueElement}
import wom.graph.CommandCallNode
import wom.graph.GraphNodePort.OutputPort
import wom.values.WomArray.WomArrayLike
import wom.values._

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

  val fullyQualifiedInputs: Map[String, WomValue] = evaluatedTaskInputs map { case (declaration, value) =>
    key.call.identifier.combine(declaration.name).fullyQualifiedName.value -> value
  }

  def findInputFilesByParameterMeta(filter: MetaValueElement => Boolean): Set[WomFile] = evaluatedTaskInputs.collect {
    case (declaration, value) if declaration.parameterMeta.exists(filter) => findFiles(value)
  }.flatten.toSet

  def findFiles(v: WomValue): Set[WomFile] = v match {
    case value: WomFile => Set(value)
    case WomOptionalValue(_, Some(value)) => findFiles(value)
    case value: WomObjectLike => value.values.values.toSet flatMap findFiles
    case WomArrayLike(value) => value.value.toSet flatMap findFiles
    case WomPair(left, right) => findFiles(left) ++ findFiles(right)
    case WomMap(_, innerMap) => (innerMap.keySet flatMap findFiles) ++ (innerMap.values.toSet flatMap findFiles)
    case _ => Set.empty
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
  
  // So it can be overridden in tests
  private [backend] lazy val cromwellFileSystems = CromwellFileSystems.instance

  lazy val configuredPathBuilderFactories: Map[String, PathBuilderFactory] = {
    cromwellFileSystems.factoriesFromConfig(backendConfig).unsafe("Failed to instantiate backend filesystem")
  }

  private lazy val configuredFactoriesWithDefault = if (configuredPathBuilderFactories.values.exists(_ == DefaultPathBuilderFactory)) {
    configuredPathBuilderFactories
  } else configuredPathBuilderFactories + DefaultPathBuilderFactory.tuple

  /**
    * Creates path builders using only the configured factories.
    */
  def pathBuilders(workflowOptions: WorkflowOptions)(implicit as: ActorSystem) = {
    PathBuilderFactory.instantiatePathBuilders(configuredPathBuilderFactories.values.toList, workflowOptions)
  }

  /**
    * Creates path builders using only the configured factories + the default factory
    */
  def pathBuildersWithDefault(workflowOptions: WorkflowOptions)(implicit as: ActorSystem) = {
    PathBuilderFactory.instantiatePathBuilders(configuredFactoriesWithDefault.values.toList, workflowOptions)
  }
}

final case class AttemptedLookupResult(name: String, value: Try[WomValue]) {
  def toPair = name -> value
}
