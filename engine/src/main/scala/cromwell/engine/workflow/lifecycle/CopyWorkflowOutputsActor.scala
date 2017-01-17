package cromwell.engine.workflow.lifecycle

import java.nio.file.Path

import akka.actor.{ActorRef, Props}
import cromwell.backend.BackendWorkflowFinalizationActor.{FinalizationFailed, FinalizationResponse, FinalizationSuccess}
import cromwell.backend.{AllBackendInitializationData, BackendConfigurationDescriptor, BackendInitializationData, BackendLifecycleActorFactory}
import cromwell.core.Dispatcher.IoDispatcher
import cromwell.core.WorkflowOptions._
import cromwell.core._
import cromwell.core.path.{PathCopier, PathFactory}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.backend.{BackendConfiguration, CromwellBackends}
import cromwell.services.io.AsyncIo
import wdl4s.values.{WdlArray, WdlMap, WdlSingleFile, WdlValue}

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps

object CopyWorkflowOutputsActor {
  def props(workflowId: WorkflowId, workflowDescriptor: EngineWorkflowDescriptor, workflowOutputs: CallOutputs,
            initializationData: AllBackendInitializationData, serviceRegistryActor: ActorRef) = Props(
    new CopyWorkflowOutputsActor(workflowId, workflowDescriptor, workflowOutputs, initializationData, serviceRegistryActor)
  ).withDispatcher(IoDispatcher)
}

class CopyWorkflowOutputsActor(workflowId: WorkflowId, val workflowDescriptor: EngineWorkflowDescriptor, workflowOutputs: CallOutputs,
                               initializationData: AllBackendInitializationData, override val serviceRegistryActor: ActorRef)
  extends EngineWorkflowFinalizationActor with PathFactory with AsyncIo {
  import cats.implicits._
  
  override val pathBuilders = workflowDescriptor.pathBuilders
  implicit val ec = context.dispatcher

  private def copyWorkflowOutputs(workflowOutputsFilePath: String): Future[Unit] = {
    val workflowOutputsPath = buildPath(workflowOutputsFilePath)

    val outputFilePaths = getOutputFilePaths

    val asyncCopies = outputFilePaths flatMap {
      case (workflowRootPath, srcPaths) =>
        // WARNING: PathCopier does not do atomic copies. The files may be partially written.
        srcPaths map { srcPath =>
          val destinationFilePath = PathCopier.getDestinationFilePath(workflowRootPath, srcPath, workflowOutputsPath)
          copy(srcPath, destinationFilePath)
        }
    }
    
    asyncCopies.toList.sequence[Future[Unit], Unit].void
  }

  private def findFiles(values: Seq[WdlValue]): Seq[WdlSingleFile] = {
    values flatMap {
      case file: WdlSingleFile => Seq(file)
      case array: WdlArray => findFiles(array.value)
      case map: WdlMap => findFiles(map.value.values.toSeq)
      case _ => Seq.empty
    }
  }
  
  private def getOutputFilePaths: Map[Path, Seq[Path]] = {
    for {
      // NOTE: Without .toSeq, outputs in arrays only yield the last output
      backend <- workflowDescriptor.backendAssignments.values.toSeq
      config <- BackendConfiguration.backendConfigurationDescriptor(backend).toOption.toSeq
      rootPath <- getBackendRootPath(backend, config).toSeq
      outputFiles = findFiles(workflowOutputs.values.map(_.wdlValue).toSeq)
      wdlPath = outputFiles map { wdlFile => rootPath.getFileSystem.getPath(wdlFile.value) }
    } yield (rootPath, wdlPath)
  } toMap

  private def getBackendRootPath(backend: String, config: BackendConfigurationDescriptor): Option[Path] = {
    getBackendFactory(backend) map getRootPath(config, initializationData.get(backend))
  }

  private def getBackendFactory(backend: String): Option[BackendLifecycleActorFactory] = {
    CromwellBackends.backendLifecycleFactoryActorByName(backend).toOption
  }

  private def getRootPath(config: BackendConfigurationDescriptor, initializationData: Option[BackendInitializationData])
                         (backendFactory: BackendLifecycleActorFactory): Path = {
    backendFactory.getExecutionRootPath(workflowDescriptor.backendDescriptor, config.backendConfig, initializationData)
  }

  final override def afterAll()(implicit ec: ExecutionContext): Future[FinalizationResponse] = {
    workflowDescriptor.getWorkflowOption(FinalWorkflowOutputsDir) map copyWorkflowOutputs getOrElse Future.successful(()) map {
      _ => FinalizationSuccess
    } recover {
      case failure => FinalizationFailed(failure)
    }
  }
}
