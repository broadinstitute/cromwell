package cromwell.engine.workflow.lifecycle

import akka.actor.{ActorRef, Props}
import cromwell.backend.BackendWorkflowFinalizationActor.{FinalizationResponse, FinalizationSuccess}
import cromwell.backend.{AllBackendInitializationData, BackendConfigurationDescriptor, BackendInitializationData, BackendLifecycleActorFactory}
import cromwell.core.Dispatcher.IoDispatcher
import cromwell.core.WorkflowOptions._
import cromwell.core._
import cromwell.core.io.AsyncIo
import cromwell.core.path.{Path, PathCopier, PathFactory}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.backend.{BackendConfiguration, CromwellBackends}
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import wdl4s.values.{WdlArray, WdlMap, WdlSingleFile, WdlValue}

import scala.concurrent.{ExecutionContext, Future}

object CopyWorkflowOutputsActor {
  def props(workflowId: WorkflowId, ioActor: ActorRef, workflowDescriptor: EngineWorkflowDescriptor, workflowOutputs: CallOutputs,
            initializationData: AllBackendInitializationData) = Props(
    new CopyWorkflowOutputsActor(workflowId, ioActor, workflowDescriptor, workflowOutputs, initializationData)
  ).withDispatcher(IoDispatcher)
}

class CopyWorkflowOutputsActor(workflowId: WorkflowId, override val ioActor: ActorRef, val workflowDescriptor: EngineWorkflowDescriptor, workflowOutputs: CallOutputs,
                               initializationData: AllBackendInitializationData)
  extends EngineWorkflowFinalizationActor with PathFactory with AsyncIo with GcsBatchCommandBuilder {

  implicit val ec = context.dispatcher
  override val pathBuilders = workflowDescriptor.pathBuilders

  override def receive = ioReceive orElse super.receive
  
  private def copyWorkflowOutputs(workflowOutputsFilePath: String): Future[Seq[Unit]] = {
    val workflowOutputsPath = buildPath(workflowOutputsFilePath)

    val outputFilePaths = getOutputFilePaths(workflowOutputsPath)

    val copies = outputFilePaths map {
      case (srcPath, dstPath) => 
        dstPath.createDirectories()
        copyAsync(srcPath, dstPath, overwrite = true)
    }
    
    Future.sequence(copies)
  }

  private def findFiles(values: Seq[WdlValue]): Seq[WdlSingleFile] = {
    values flatMap {
      case file: WdlSingleFile => Seq(file)
      case array: WdlArray => findFiles(array.value)
      case map: WdlMap => findFiles(map.value.values.toSeq)
      case _ => Seq.empty
    }
  }
  
  private def getOutputFilePaths(workflowOutputsPath: Path): List[(Path, Path)] = {
    val rootAndFiles = for {
      // NOTE: Without .toSeq, outputs in arrays only yield the last output
      backend <- workflowDescriptor.backendAssignments.values.toSeq
      config <- BackendConfiguration.backendConfigurationDescriptor(backend).toOption.toSeq
      rootPath <- getBackendRootPath(backend, config).toSeq
      outputFiles = findFiles(workflowOutputs.values.map(_.wdlValue).toSeq).map(_.value)
    } yield (rootPath, outputFiles)
    
    val outputFileDestinations = rootAndFiles flatMap {
      case (workflowRoot, outputs) =>
        outputs map { output => 
          val outputPath = PathFactory.buildPath(output, pathBuilders)
          outputPath -> PathCopier.getDestinationFilePath(workflowRoot, outputPath, workflowOutputsPath) 
        }
    }
    outputFileDestinations.distinct.toList
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

  final override def afterAll()(implicit ec: ExecutionContext): Future[FinalizationResponse] = {
    workflowDescriptor.getWorkflowOption(FinalWorkflowOutputsDir) match {
      case Some(outputs) => copyWorkflowOutputs(outputs) map { _ => FinalizationSuccess }
      case None => Future.successful(FinalizationSuccess)
    }
  }
}
