package cromwell.backend.impl.jes

import java.nio.file.Path

import akka.actor.Props
import better.files._
import cromwell.backend.impl.jes.io._
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor, BackendWorkflowFinalizationActor}
import cromwell.core.{ExecutionStore, OutputStore, PathCopier}
import wdl4s.Call

import scala.concurrent.Future

object JesFinalizationActor {
  def props(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], jesConfiguration: JesConfiguration,
            executionStore: ExecutionStore, outputStore: OutputStore) = {
    Props(new JesFinalizationActor(workflowDescriptor, calls, jesConfiguration, executionStore, outputStore))
      .withDispatcher("akka.dispatchers.slow-actor-dispatcher")
  }
}

class JesFinalizationActor (override val workflowDescriptor: BackendWorkflowDescriptor,
                            override val calls: Seq[Call],
                            jesConfiguration: JesConfiguration, executionStore: ExecutionStore,
                            outputStore: OutputStore) extends BackendWorkflowFinalizationActor {

  override val configurationDescriptor = jesConfiguration.configurationDescriptor

  private val workflowPaths = new JesWorkflowPaths(workflowDescriptor, jesConfiguration)

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
    import cromwell.core.WorkflowOptions._
    /*
    NOTE: Only using one thread pool slot here to upload all the files for all the calls.
    Using the slow-actor-dispatcher defined in application.conf because this might take a while.
    One could also use Future.sequence to flood the dispatcher, or even create a separate jes final call specific thread
    pool for parallel uploads.

    Measure and optimize as necessary. Will likely need retry code at some level as well.
     */
    Future(workflowDescriptor.getWorkflowOption(FinalCallLogsDir) foreach copyCallOutputs)(context.system.dispatcher)
  }

  private def copyCallOutputs(callLogsDir: String): Unit = {
    copyLogs(toJesPath(callLogsDir), logPaths)
  }

  private def toJesPath(value: String): Path = {
    val fileSystem = buildFilesystem(workflowDescriptor,
      jesConfiguration.jesAttributes.gcsFilesystemAuth, jesConfiguration.googleConfig)

    fileSystem.getPath(value)
  }

  private lazy val logPaths: Seq[Path] = {
    val allCallPaths = executionStore.store.toSeq collect {
      case (backendJobDescriptorKey: BackendJobDescriptorKey, _) =>
        JesCallPaths(backendJobDescriptorKey, workflowDescriptor, jesConfiguration)
    }

    allCallPaths flatMap { callPaths =>
      Seq(callPaths.stdoutPath, callPaths.stderrPath, callPaths.jesLogPath)
    }
  }

  private def getWorkflowOption(key: String): Option[String] = {
    val workflowOptions = workflowDescriptor.workflowOptions
    workflowOptions.get(key).toOption
  }

  private def copyLogs(callLogsDirPath: Path, logPaths: Seq[Path]): Unit = {
    logPaths.foreach(PathCopier.copy(workflowPaths.rootPath, _, callLogsDirPath))
  }
}
