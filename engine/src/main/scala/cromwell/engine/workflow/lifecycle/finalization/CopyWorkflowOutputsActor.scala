package cromwell.engine.workflow.lifecycle.finalization

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import akka.event.LoggingReceive
import cromwell.backend.BackendLifecycleActor.BackendWorkflowLifecycleActorResponse
import cromwell.backend.BackendWorkflowFinalizationActor.{FinalizationFailed, FinalizationResponse, FinalizationSuccess, Finalize}
import cromwell.backend.{AllBackendInitializationData, BackendConfigurationDescriptor, BackendInitializationData, BackendLifecycleActorFactory}
import cromwell.core.Dispatcher.IoDispatcher
import cromwell.core.WorkflowOptions._
import cromwell.core._
import cromwell.core.io.AsyncIoActorClient
import cromwell.core.path.{Path, PathCopier, PathFactory}
import cromwell.engine.EngineWorkflowDescriptor
import cromwell.engine.backend.{BackendConfiguration, CromwellBackends}
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import wom.values.{WomSingleFile, WomValue}

import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Success}

object CopyWorkflowOutputsActor {
  def props(workflowId: WorkflowId, ioActor: ActorRef, workflowDescriptor: EngineWorkflowDescriptor, workflowOutputs: CallOutputs,
            initializationData: AllBackendInitializationData) = Props(
    new CopyWorkflowOutputsActor(workflowId, ioActor, workflowDescriptor, workflowOutputs, initializationData)
  ).withDispatcher(IoDispatcher)
}

class CopyWorkflowOutputsActor(workflowId: WorkflowId, override val ioActor: ActorRef, val workflowDescriptor: EngineWorkflowDescriptor, workflowOutputs: CallOutputs,
                               initializationData: AllBackendInitializationData)
  extends Actor with ActorLogging with PathFactory with AsyncIoActorClient {
  override lazy val ioCommandBuilder = GcsBatchCommandBuilder
  implicit val ec = context.dispatcher
  override val pathBuilders = workflowDescriptor.pathBuilders

  override def receive = LoggingReceive {
    case Finalize => performActionThenRespond(afterAll()(context.dispatcher), FinalizationFailed)(context.dispatcher)
  }

  private def performActionThenRespond(operation: => Future[BackendWorkflowLifecycleActorResponse],
                                       onFailure: (Throwable) => BackendWorkflowLifecycleActorResponse)
                                      (implicit ec: ExecutionContext) = {
    val respondTo: ActorRef = sender
    operation onComplete {
      case Success(r) => respondTo ! r
      case Failure(t) => respondTo ! onFailure(t)
    }
  }

  private def copyWorkflowOutputs(workflowOutputsFilePath: String): Future[Seq[Unit]] = {
    val workflowOutputsPath = buildPath(workflowOutputsFilePath)

    val outputFilePaths = getOutputFilePaths(workflowOutputsPath)

    val copies = outputFilePaths map {
      case (srcPath, dstPath) => 
        dstPath.createDirectories()
        asyncIo.copyAsync(srcPath, dstPath)
    }
    
    Future.sequence(copies)
  }

  private def findFiles(values: Seq[WomValue]): Seq[WomSingleFile] = {
    values flatMap {
      _.collectAsSeq {
        case file: WomSingleFile => file
      }
    }
  }
  
  private def getOutputFilePaths(workflowOutputsPath: Path): List[(Path, Path)] = {
    val rootAndFiles = for {
      // NOTE: Without .toSeq, outputs in arrays only yield the last output
      backend <- workflowDescriptor.backendAssignments.values.toSeq
      config <- BackendConfiguration.backendConfigurationDescriptor(backend).toOption.toSeq
      rootPath <- getBackendRootPath(backend, config).toSeq
      outputFiles = findFiles(workflowOutputs.outputs.values.toSeq).map(_.value)
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

  /**
    * Happens after everything else runs
    */
  final def afterAll()(implicit ec: ExecutionContext): Future[FinalizationResponse] = {
    workflowDescriptor.getWorkflowOption(FinalWorkflowOutputsDir) match {
      case Some(outputs) => copyWorkflowOutputs(outputs) map { _ => FinalizationSuccess }
      case None => Future.successful(FinalizationSuccess)
    }
  }
}
