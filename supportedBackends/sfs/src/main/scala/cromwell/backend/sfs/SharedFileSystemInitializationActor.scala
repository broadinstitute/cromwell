package cromwell.backend.sfs

import akka.actor.ActorRef
import better.files._
import cromwell.backend.io.{WorkflowPaths, WorkflowPathsBackendInitializationData}
import cromwell.backend.wfs.{DefaultWorkflowFileSystemProvider, WorkflowFileSystemProvider}
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor, BackendWorkflowInitializationActor}
import cromwell.core.Dispatcher
import wdl4s.{Call, WdlExpression}

import scala.concurrent.Future
import scala.util.Try

case class SharedFileSystemInitializationActorParams
(
  serviceRegistryActor: ActorRef,
  workflowDescriptor: BackendWorkflowDescriptor,
  configurationDescriptor: BackendConfigurationDescriptor,
  calls: Seq[Call],
  supportsDocker: Boolean,
  runtimeAttributesBuilder: SharedFileSystemValidatedRuntimeAttributesBuilder
)

case class SharedFileSystemBackendInitializationData(workflowPaths: WorkflowPaths)
  extends WorkflowPathsBackendInitializationData

class SharedFileSystemInitializationActor(params: SharedFileSystemInitializationActorParams)
  extends BackendWorkflowInitializationActor {

  override lazy val workflowDescriptor: BackendWorkflowDescriptor = params.workflowDescriptor
  override lazy val configurationDescriptor: BackendConfigurationDescriptor = params.configurationDescriptor
  override lazy val calls: Seq[Call] = params.calls
  override lazy val serviceRegistryActor: ActorRef = params.serviceRegistryActor

  lazy val runtimeAttributesBuilder: SharedFileSystemValidatedRuntimeAttributesBuilder =
    params.runtimeAttributesBuilder.withDockerSupport(params.supportsDocker)

  override protected def runtimeAttributeValidators: Map[String, (Option[WdlExpression]) => Boolean] = {
    runtimeAttributesBuilder.validations.map(validation =>
      validation.key -> validation.validateOptionalExpression _
    ).toMap
  }

  val providers = Seq(GcsWorkflowFileSystemProvider, DefaultWorkflowFileSystemProvider)
  val ioDispatcher = context.system.dispatchers.lookup(Dispatcher.IoDispatcher)

  val workflowPaths = WorkflowFileSystemProvider.workflowPaths(configurationDescriptor, workflowDescriptor,
    providers, ioDispatcher)

  /**
    * A call which happens before anything else runs
    */
  override def beforeAll(): Future[Option[BackendInitializationData]] = {
    Future.fromTry(Try {
      publishWorkflowRoot(workflowPaths.workflowRoot.toString)
      workflowPaths.workflowRoot.createDirectories()
      Option(SharedFileSystemBackendInitializationData(workflowPaths))
    })
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
}
