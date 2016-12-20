package cromwell.backend.impl.jes

import java.nio.file.Path

import akka.actor.ActorRef
import cromwell.backend.standard.{StandardFinalizationActor, StandardFinalizationActorParams}
import cromwell.backend.{BackendWorkflowDescriptor, JobExecutionMap, _}
import cromwell.core.CallOutputs
import cromwell.core.path.PathCopier
import cromwell.services.io.AsyncIo
import wdl4s.TaskCall

import scala.concurrent.Future

case class JesFinalizationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  calls: Set[TaskCall],
  serviceRegistryActor: ActorRef,
  jesConfiguration: JesConfiguration,
  jobExecutionMap: JobExecutionMap,
  workflowOutputs: CallOutputs,
  initializationDataOption: Option[BackendInitializationData]
) extends StandardFinalizationActorParams {
  override val configurationDescriptor: BackendConfigurationDescriptor = jesConfiguration.configurationDescriptor
}

class JesFinalizationActor(jesParams: JesFinalizationActorParams)
  extends StandardFinalizationActor with AsyncIo {

  override val standardParams: StandardFinalizationActorParams = jesParams

  lazy val jesConfiguration: JesConfiguration = jesParams.jesConfiguration

  private val workflowPaths = initializationDataOption.map {
    _.asInstanceOf[JesBackendInitializationData].workflowPaths
  }

  override def afterAll(): Future[Unit] = {
    for {
      // NOTE: These are currently in series, not in parallel. Not sure how many threads to throw at finalization
      _ <- deleteAuthenticationFile()
      _ <- copyCallOutputs()
    } yield ()
  }

  private def deleteAuthenticationFile(): Future[Unit] = {
    (jesConfiguration.needAuthFileUpload, workflowPaths) match {
      case (true, Some(paths)) => delete(paths.gcsAuthFilePath, swallowIOExceptions = false)
      case _ => Future.successful(())
    }
  }

  private def copyCallOutputs(): Future[Unit] = {
    /*
    NOTE: All copies are sent in parallel to the ioActor. 
    We still wait for all of them to be complete before returning the Future (see Future.sequence)
    We could make this even more async and have an actor waiting for all the futures to come back.
    Or even better create a dedicated actor and use the ack based version of the Io messages.
     */
    workflowPaths match {
      case Some(paths) => paths.finalCallLogsPath map copyCallOutputs getOrElse Future.successful(())
      case _ => Future.successful(())
    }
  }

  private def copyCallOutputs(callLogsPath: Path): Future[Unit] = {
    copyLogs(callLogsPath, logPaths)
  }

  private lazy val logPaths: Seq[Path] = {
    val allCallPaths = jobExecutionMap flatMap {
      case (backendJobDescriptor, keys) =>
        keys map { JesJobPaths(_, backendJobDescriptor, jesConfiguration)(context.system) }
    }

    allCallPaths.toSeq flatMap { callPaths =>
      Seq(callPaths.stdout, callPaths.stderr, callPaths.jesLogPath)
    }
  }

  private def copyLogs(callLogsDirPath: Path, logPaths: Seq[Path]): Future[Unit] = {
    workflowPaths match {
      case Some(paths) => 
        val logAsyncCopies = logPaths map { sourceFilePath =>
          val destinationPath = PathCopier.getDestinationFilePath(paths.executionRoot, sourceFilePath, callLogsDirPath)
          copy(sourceFilePath, destinationPath)
        }
        Future.sequence(logAsyncCopies) map { _ => () }
      case None => Future.successful(())
    }
  }
}
