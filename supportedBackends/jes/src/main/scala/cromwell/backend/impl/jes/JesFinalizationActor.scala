package cromwell.backend.impl.jes

import java.nio.file.Path

import better.files._
import cats.instances.future._
import cats.syntax.functor._
import cromwell.backend._
import cromwell.backend.standard.{StandardFinalizationActor, StandardFinalizationActorParams}
import cromwell.core.CallOutputs
import cromwell.core.Dispatcher.IoDispatcher
import cromwell.core.path.PathCopier
import wdl4s.TaskCall

import scala.concurrent.Future
import scala.language.postfixOps

case class JesFinalizationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  calls: Set[TaskCall],
  jesConfiguration: JesConfiguration,
  jobExecutionMap: JobExecutionMap,
  workflowOutputs: CallOutputs,
  initializationDataOption: Option[BackendInitializationData]
) extends StandardFinalizationActorParams {
  override val configurationDescriptor: BackendConfigurationDescriptor = jesConfiguration.configurationDescriptor
}

class JesFinalizationActor(jesParams: JesFinalizationActorParams)
  extends StandardFinalizationActor {

  override val standardParams: StandardFinalizationActorParams = jesParams

  lazy val jesConfiguration: JesConfiguration = jesParams.jesConfiguration

  private val workflowPaths = initializationDataOption.map {
    _.asInstanceOf[JesBackendInitializationData].workflowPaths
  }

  private val iOExecutionContext = context.system.dispatchers.lookup(IoDispatcher)

  override def afterAll(): Future[Unit] = {
    for {
      // NOTE: These are currently in series, not in parallel. Not sure how many threads to throw at finalization
      _ <- deleteAuthenticationFile()
      _ <- copyCallOutputs()
    } yield ()
  }

  private def deleteAuthenticationFile(): Future[Unit] = {
    (jesConfiguration.needAuthFileUpload, workflowPaths) match {
      case (true, Some(paths)) => Future { File(paths.gcsAuthFilePath).delete(false) } void
      case _ => Future.successful(())
    }
  }

  private def copyCallOutputs(): Future[Unit] = {
    /*
    NOTE: Only using one thread pool slot here to upload all the files for all the calls.
    Using the io-dispatcher defined in application.conf because this might take a while.
    One could also use Future.sequence to flood the dispatcher, or even create a separate jes final call specific thread
    pool for parallel uploads.

    Measure and optimize as necessary. Will likely need retry code at some level as well.
     */
    workflowPaths match {
      case Some(paths) => Future(paths.finalCallLogsPath foreach copyCallOutputs)(iOExecutionContext)
      case _ => Future.successful(())
    }
  }

  private def copyCallOutputs(callLogsPath: Path): Unit = {
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

  private def copyLogs(callLogsDirPath: Path, logPaths: Seq[Path]): Unit = {
    workflowPaths match {
      case Some(paths) => logPaths.foreach(PathCopier.copy(paths.executionRoot, _, callLogsDirPath))
      case None =>
    }
  }
}
