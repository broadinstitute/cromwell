package cromwell.engine.workflow.lifecycle

import java.nio.file.Path

import akka.actor.Props
import cromwell.backend.BackendWorkflowFinalizationActor.{FinalizationResponse, FinalizationSuccess}
import cromwell.backend.{AllBackendInitializationData, BackendConfigurationDescriptor, BackendInitializationData, BackendLifecycleActorFactory}
import cromwell.core.Dispatcher.IoDispatcher
import cromwell.core.WorkflowOptions._
import cromwell.core._
import cromwell.core.path.{PathCopier, PathFactory}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.backend.{BackendConfiguration, CromwellBackends}
import wdl4s.values.{WdlArray, WdlMap, WdlSingleFile, WdlValue}

import scala.concurrent.{ExecutionContext, Future}

object CopyWorkflowOutputsActor {
  def props(workflowId: WorkflowId, workflowDescriptor: EngineWorkflowDescriptor, workflowOutputs: CallOutputs,
            initializationData: AllBackendInitializationData) = Props(
    new CopyWorkflowOutputsActor(workflowId, workflowDescriptor, workflowOutputs, initializationData)
  ).withDispatcher(IoDispatcher)
}

class CopyWorkflowOutputsActor(workflowId: WorkflowId, val workflowDescriptor: EngineWorkflowDescriptor, workflowOutputs: CallOutputs,
                               initializationData: AllBackendInitializationData)
  extends EngineWorkflowFinalizationActor with PathFactory {

  override val pathBuilders = workflowDescriptor.pathBuilders

  private def copyWorkflowOutputs(workflowOutputsFilePath: String): Unit = {
    val workflowOutputsPath = buildPath(workflowOutputsFilePath)

    val outputFilePaths = getOutputFilePaths

    outputFilePaths foreach {
      case (workflowRootPath, srcPath) =>
        // WARNING: PathCopier does not do atomic copies. The files may be partially written.
        PathCopier.copy(workflowRootPath, srcPath, workflowOutputsPath)
    }
  }

  private def findFiles(values: Seq[WdlValue]): Seq[WdlSingleFile] = {
    values flatMap {
      case file: WdlSingleFile => Seq(file)
      case array: WdlArray => findFiles(array.value)
      case map: WdlMap => findFiles(map.value.values.toSeq)
      case _ => Seq.empty
    }
  }
  
  private def getOutputFilePaths: Seq[(Path, Path)] = {
    for {
      // NOTE: Without .toSeq, outputs in arrays only yield the last output
      backend <- workflowDescriptor.backendAssignments.values.toSeq
      config <- BackendConfiguration.backendConfigurationDescriptor(backend).toOption.toSeq
      rootPath <- getBackendRootPath(backend, config).toSeq
      outputFiles = findFiles(workflowOutputs.values.map(_.wdlValue).toSeq)
      wdlFile <- outputFiles
      wdlPath = rootPath.getFileSystem.getPath(wdlFile.value)
    } yield (rootPath, wdlPath)
  }

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

  final override def afterAll()(implicit ec: ExecutionContext): Future[FinalizationResponse] = Future {
    workflowDescriptor.getWorkflowOption(FinalWorkflowOutputsDir) foreach copyWorkflowOutputs
    FinalizationSuccess
  }
}
