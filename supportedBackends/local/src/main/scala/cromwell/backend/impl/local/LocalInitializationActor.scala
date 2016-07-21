package cromwell.backend.impl.local

import java.nio.file.FileSystems

import akka.actor.{ActorRef, Props}
import better.files._
import cromwell.backend.impl.local.LocalInitializationActor._
import cromwell.backend.io.WorkflowPaths
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor, BackendWorkflowInitializationActor}
import cromwell.filesystems.gcs.{GcsFileSystemProvider, GcsFileSystem}
import wdl4s.types.{WdlBooleanType, WdlStringType}
import wdl4s.{Call, WdlExpression}
import LocalImplicits._

import scala.concurrent.Future

object LocalInitializationActor {
  val SupportedKeys = Set(DockerKey, FailOnStderrKey, ContinueOnReturnCodeKey)

  def props(workflowDescriptor: BackendWorkflowDescriptor,
            calls: Seq[Call],
            configurationDescriptor: BackendConfigurationDescriptor,
            serviceRegistryActor: ActorRef,
            localConfiguration: LocalConfiguration): Props =
    Props(new LocalInitializationActor(workflowDescriptor, calls, configurationDescriptor, serviceRegistryActor, localConfiguration))

}

class LocalInitializationActor(override val workflowDescriptor: BackendWorkflowDescriptor,
                               override val calls: Seq[Call],
                               override val configurationDescriptor: BackendConfigurationDescriptor,
                               override val serviceRegistryActor: ActorRef,
                               localConfiguration: LocalConfiguration) extends BackendWorkflowInitializationActor {

  override protected def runtimeAttributeValidators: Map[String, (Option[WdlExpression]) => Boolean] = Map(
    DockerKey -> wdlTypePredicate(valueRequired = false, WdlStringType.isCoerceableFrom),
    FailOnStderrKey -> wdlTypePredicate(valueRequired = false, WdlBooleanType.isCoerceableFrom),
    ContinueOnReturnCodeKey -> continueOnReturnCodePredicate(valueRequired = false)
  )

  private val fileSystems = {
    val maybeGcs = {
      localConfiguration.gcsAuthMode map { authMode =>
        val storage = authMode.buildStorage(workflowDescriptor.workflowOptions.toGoogleAuthOptions, localConfiguration.googleConfig)
        GcsFileSystem(GcsFileSystemProvider(storage))
      }
    }

    List(maybeGcs, Option(FileSystems.getDefault)).flatten
  }

  private val workflowPaths = new WorkflowPaths(workflowDescriptor, configurationDescriptor.backendConfig, fileSystems)

  /**
    * A call which happens before anything else runs
    */
  override def beforeAll(): Future[Option[BackendInitializationData]] = {
    publishWorkflowRoot(workflowPaths.workflowRoot.toString)
    workflowPaths.workflowRoot.createDirectories()
    Future.successful(Option(LocalBackendInitializationData(workflowPaths)))
  }

  /**
    * Log a warning if there are non-supported runtime attributes defined for the call.
    */
  override def validate(): Future[Unit] = {
    Future {
      calls foreach { call =>
        val runtimeAttributes = call.task.runtimeAttributes.attrs
        val notSupportedAttributes = runtimeAttributes filterKeys { !SupportedKeys.contains(_) }

        if (notSupportedAttributes.nonEmpty) {
          val notSupportedAttrString = notSupportedAttributes.keys mkString ", "
          workflowLogger.warn(s"Key/s [$notSupportedAttrString] is/are not supported by LocalBackend. Unsupported attributes will not be part of jobs executions.")
        }
      }
    }
  }
}
