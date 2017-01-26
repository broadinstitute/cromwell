package cromwell.backend.impl.tes

import akka.actor.{ActorRef, Props}
import better.files.File
import cromwell.backend.standard._
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendWorkflowDescriptor}
import cromwell.core.Dispatcher.BackendDispatcher
import cromwell.core.path.FileImplicits._
import cromwell.core.path.{DefaultPathBuilderFactory, PathBuilder, PathBuilderFactory}
import wdl4s.TaskCall

import scala.concurrent.Future
import scala.util.Try

object TesInitializationActor {
  /* NOTE: Only used by tests */
  def props(workflowDescriptor: BackendWorkflowDescriptor,
            calls: Set[TaskCall],
            tesConfiguration: TesConfiguration,
            serviceRegistryActor: ActorRef): Props = {
    val params = TesInitializationActorParams(workflowDescriptor, calls, tesConfiguration, serviceRegistryActor)
    Props(new TesInitializationActor(params)).withDispatcher(BackendDispatcher)
  }
}

case class TesInitializationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  calls: Set[TaskCall],
  tesConfiguration: TesConfiguration,
  serviceRegistryActor: ActorRef
) extends StandardInitializationActorParams {
  override val configurationDescriptor: BackendConfigurationDescriptor = tesConfiguration.configurationDescriptor
}

class TesInitializationActor(params: TesInitializationActorParams)
  extends StandardInitializationActor {

  override val standardParams: StandardInitializationActorParams = params

  private val tesConfiguration = params.tesConfiguration

  lazy val pathBuilderFactories: List[PathBuilderFactory] = List(Option(DefaultPathBuilderFactory)).flatten

  lazy val pathBuilders: List[PathBuilder] =
    pathBuilderFactories map {
      _.withOptions(workflowDescriptor.workflowOptions)(context.system)
    }

  val workflowPaths: TesWorkflowPaths =
    new TesWorkflowPaths(workflowDescriptor, tesConfiguration.configurationDescriptor.backendConfig, pathBuilders)

  override lazy val runtimeAttributesBuilder: StandardValidatedRuntimeAttributesBuilder =
    TesRuntimeAttributes.runtimeAttributesBuilder

  override def beforeAll(): Future[Option[BackendInitializationData]] = {
    Future.fromTry(Try {
      publishWorkflowRoot(workflowPaths.workflowRoot.toString)
      File(workflowPaths.workflowRoot).createPermissionedDirectories()
      Option(TesBackendInitializationData(workflowPaths, runtimeAttributesBuilder, tesConfiguration))
    })
  }
}
