package cromwell.backend.impl.jes

import java.nio.file.Path

import akka.actor.Props
import better.files._
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor, BackendWorkflowFinalizationActor}
import cromwell.core.{ExecutionStore, OutputStore, PathCopier}
import wdl4s.Call

import scala.concurrent.Future

object JesFinalizationActor {
  def props(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], jesConfiguration: JesConfiguration,
            executionStore: ExecutionStore, outputStore: OutputStore, initializationData: JesBackendInitializationData) = {
    Props(new JesFinalizationActor(workflowDescriptor, calls, jesConfiguration, executionStore, outputStore, initializationData))
  }
}

class JesFinalizationActor (override val workflowDescriptor: BackendWorkflowDescriptor,
                            override val calls: Seq[Call],
                            jesConfiguration: JesConfiguration, executionStore: ExecutionStore,
                            outputStore: OutputStore,
                            initializationData: JesBackendInitializationData) extends BackendWorkflowFinalizationActor {

  override val configurationDescriptor = jesConfiguration.configurationDescriptor

  private val workflowPaths = initializationData.workflowPaths

  private val iOExecutionContext = context.system.dispatchers.lookup("akka.dispatchers.io-dispatcher")

  override def afterAll(): Future[Unit] = {
    for {
      // NOTE: These are currently in series, not in parallel. Not sure how many threads to throw at finalization
      _ <- deleteAuthenticationFile()
      _ <- copyCallOutputs()
    } yield ()
  }

  private def deleteAuthenticationFile(): Future[Unit] = {
    if (jesConfiguration.needAuthFileUpload) {
      Future(workflowPaths.gcsAuthFilePath.delete(false)) map { _ => () }
    } else {
      Future.successful(())
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
    Future(workflowPaths.finalCallLogsPath foreach copyCallOutputs)(iOExecutionContext)
  }

  private def copyCallOutputs(callLogsPath: Path): Unit = {
    copyLogs(callLogsPath, logPaths)
  }

  private lazy val logPaths: Seq[Path] = {
    val allCallPaths = executionStore.store.toSeq collect {
      case (backendJobDescriptorKey: BackendJobDescriptorKey, _) =>
        initializationData.workflowPaths.toJesCallPaths(backendJobDescriptorKey)
    }

    allCallPaths flatMap { callPaths =>
      Seq(callPaths.stdoutPath, callPaths.stderrPath, callPaths.jesLogPath)
    }
  }

  private def copyLogs(callLogsDirPath: Path, logPaths: Seq[Path]): Unit = {
    logPaths.foreach(PathCopier.copy(workflowPaths.rootPath, _, callLogsDirPath))
  }
}
