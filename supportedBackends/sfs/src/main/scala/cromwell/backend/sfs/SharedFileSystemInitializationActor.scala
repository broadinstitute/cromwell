package cromwell.backend.sfs

import akka.actor.ActorRef
import better.files._
import cromwell.backend.io.{WorkflowPaths, WorkflowPathsBackendInitializationData}
import cromwell.backend.validation.RuntimeAttributesDefault
import cromwell.backend.wfs.WorkflowPathBuilder
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor, BackendWorkflowInitializationActor}
import cromwell.core.WorkflowOptions
import cromwell.core.path.PathBuilderFactory
import wdl4s.TaskCall
import wdl4s.values.WdlValue

import scala.concurrent.Future
import scala.util.Try

case class SharedFileSystemInitializationActorParams
(
  serviceRegistryActor: ActorRef,
  workflowDescriptor: BackendWorkflowDescriptor,
  configurationDescriptor: BackendConfigurationDescriptor,
  calls: Set[TaskCall],
  pathBuilderFactories: List[PathBuilderFactory]
)

class SharedFileSystemBackendInitializationData
(
  val workflowPaths: WorkflowPaths,
  val runtimeAttributesBuilder: SharedFileSystemValidatedRuntimeAttributesBuilder)
  extends WorkflowPathsBackendInitializationData

/**
  * Initializes a shared file system actor factory and creates initialization data to pass to the execution actors.
  *
  * @param params Initialization parameters.
  */
class SharedFileSystemInitializationActor(params: SharedFileSystemInitializationActorParams)
  extends BackendWorkflowInitializationActor {

  override lazy val workflowDescriptor: BackendWorkflowDescriptor = params.workflowDescriptor
  override lazy val configurationDescriptor: BackendConfigurationDescriptor = params.configurationDescriptor
  override lazy val calls: Set[TaskCall] = params.calls
  override lazy val serviceRegistryActor: ActorRef = params.serviceRegistryActor

  def runtimeAttributesBuilder: SharedFileSystemValidatedRuntimeAttributesBuilder =
    SharedFileSystemValidatedRuntimeAttributesBuilder.default

  override protected def runtimeAttributeValidators: Map[String, (Option[WdlValue]) => Boolean] = {
    runtimeAttributesBuilder.validations.map(validation =>
      validation.key -> validation.validateOptionalExpression _
    ).toMap
  }

  val pathBuilders = params.pathBuilderFactories map { _.withOptions(workflowDescriptor.workflowOptions)(context.system) }

  val workflowPaths = WorkflowPathBuilder.workflowPaths(configurationDescriptor, workflowDescriptor, pathBuilders)

  override def beforeAll(): Future[Option[BackendInitializationData]] = {
    Future.fromTry(Try {
      publishWorkflowRoot(workflowPaths.workflowRoot.toString)
      File(workflowPaths.workflowRoot).createDirectories()
      Option(initializationData)
    })
  }

  def initializationData: SharedFileSystemBackendInitializationData = {
    new SharedFileSystemBackendInitializationData(workflowPaths, runtimeAttributesBuilder)
  }

  /**
    * Log a warning if there are non-supported runtime attributes defined for the call.
    */
  override def validate(): Future[Unit] = {
    Future.fromTry(Try {
      calls foreach { call =>
        val runtimeAttributeKeys = call.task.runtimeAttributes.attrs.keys.toList
        val notSupportedAttributes = runtimeAttributesBuilder.unsupportedKeys(runtimeAttributeKeys).toList

        if (notSupportedAttributes.nonEmpty) {
          val notSupportedAttrString = notSupportedAttributes mkString ", "
          workflowLogger.warn(
            s"Key/s [$notSupportedAttrString] is/are not supported by backend. " +
              s"Unsupported attributes will not be part of jobs executions.")
        }
      }
    })
  }

  override protected def coerceDefaultRuntimeAttributes(options: WorkflowOptions): Try[Map[String, WdlValue]] = {
    RuntimeAttributesDefault.workflowOptionsDefault(options, runtimeAttributesBuilder.validations.map(v => v.key -> v.coercion).toMap)
  }
}
