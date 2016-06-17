package cromwell.engine.workflow.lifecycle

import java.nio.file.Path

import akka.actor.Props
import cromwell.backend.BackendWorkflowFinalizationActor.{FinalizationSuccess, FinalizationResponse}
import cromwell.backend.{BackendConfigurationDescriptor, BackendLifecycleActorFactory}
import cromwell.core._
import cromwell.core.WorkflowOptions._
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.backend.{BackendConfiguration, CromwellBackends}
import wdl4s.ReportableSymbol
import wdl4s.values.WdlSingleFile

import scala.concurrent.{Future, ExecutionContext}

object CopyWorkflowOutputsActor {
  def props(workflowId: WorkflowId, workflowDescriptor: EngineWorkflowDescriptor, outputStore: OutputStore) = Props(
    new CopyWorkflowOutputsActor(workflowId, workflowDescriptor, outputStore)
  )
}

class CopyWorkflowOutputsActor(workflowId: WorkflowId, val workflowDescriptor: EngineWorkflowDescriptor, outputStore: OutputStore)
  extends EngineWorkflowFinalizationActor with PathFactory {

  private def copyWorkflowOutputs(workflowOutputsFilePath: String): Unit = {
    val workflowOutputsPath = buildPath(workflowOutputsFilePath, workflowDescriptor.engineFilesystems)

    val reportableOutputs = workflowDescriptor.backendDescriptor.workflowNamespace.workflow.outputs

    val outputFilePaths = getOutputFilePaths(reportableOutputs)

    outputFilePaths foreach {
      case (workflowRootPath, srcPath) =>
        PathCopier.copy(workflowRootPath, srcPath, workflowOutputsPath)
    }
  }

  private def getOutputFilePaths(reportableOutputs: Seq[ReportableSymbol]): Seq[(Path, Path)] = {
    for {
      reportableOutput <- reportableOutputs
      // NOTE: Without .toSeq, outputs in arrays only yield the last output
      (backend, calls) <- workflowDescriptor.backendAssignments.groupBy(_._2).mapValues(_.keys.toSeq).toSeq
      config <- BackendConfiguration.backendConfigurationDescriptor(backend).toOption.toSeq
      rootPath <- getBackendRootPath(backend, config).toSeq
      call <- calls
      // NOTE: Without .toSeq, outputs in arrays only yield the last output
      (outputCallKey, outputEntries) <- outputStore.store.toSeq
      // Only get paths for the original scatter call, not the indexed entries
      if outputCallKey.call == call && outputCallKey.index.isEmpty
      outputEntry <- outputEntries
      if reportableOutput.fullyQualifiedName == s"${call.fullyQualifiedName}.${outputEntry.name}"
      wdlValue <- outputEntry.wdlValue.toSeq
      collected = wdlValue collectAsSeq { case f: WdlSingleFile => f }
      wdlFile <- collected
      wdlPath = rootPath.getFileSystem.getPath(wdlFile.value)
    } yield (rootPath, wdlPath)
  }

  private def getBackendRootPath(backend: String, config: BackendConfigurationDescriptor): Option[Path] = {
    getBackendFactory(backend) map getRootPath(config)
  }

  private def getBackendFactory(backend: String): Option[BackendLifecycleActorFactory] = {
    CromwellBackends.shadowBackendLifecycleFactory(backend).toOption
  }

  private def getRootPath(config: BackendConfigurationDescriptor)
                         (backendFactory: BackendLifecycleActorFactory): Path = {
    backendFactory.getExecutionRootPath(workflowDescriptor.backendDescriptor, config.backendConfig)
  }

  final override def afterAll()(implicit ec: ExecutionContext): Future[FinalizationResponse] = Future {
    workflowDescriptor.getWorkflowOption(FinalWorkflowOutputsDir) foreach copyWorkflowOutputs
    FinalizationSuccess
  }
}
