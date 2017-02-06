package cromwell.backend.standard

import java.nio.file.Path

import akka.dispatch.MessageDispatcher
import cromwell.backend._
import cromwell.backend.io.WorkflowPaths
import cromwell.core.CallOutputs
import cromwell.core.Dispatcher.IoDispatcher
import cromwell.core.path.PathCopier
import wdl4s.TaskCall

import scala.concurrent.Future

trait StandardFinalizationActorParams {
  def workflowDescriptor: BackendWorkflowDescriptor

  def calls: Set[TaskCall]

  def jobExecutionMap: JobExecutionMap

  def workflowOutputs: CallOutputs

  def initializationDataOption: Option[BackendInitializationData]

  def configurationDescriptor: BackendConfigurationDescriptor
}

case class DefaultStandardFinalizationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  calls: Set[TaskCall],
  jobExecutionMap: JobExecutionMap,
  workflowOutputs: CallOutputs,
  initializationDataOption: Option[BackendInitializationData],
  configurationDescriptor: BackendConfigurationDescriptor
) extends StandardFinalizationActorParams

/**
  * Implements BackendWorkflowFinalizationActor.afterAll() by calling copyCallLogs.
  *
  * Subclasses should either call super.afterAll(), or invoke copyCallLogs().
  *
  * @param standardParams Standard parameters.
  */
class StandardFinalizationActor(val standardParams: StandardFinalizationActorParams)
  extends BackendWorkflowFinalizationActor {

  override lazy val workflowDescriptor: BackendWorkflowDescriptor = standardParams.workflowDescriptor
  override lazy val calls: Set[TaskCall] = standardParams.calls
  lazy val initializationDataOption: Option[BackendInitializationData] = standardParams.initializationDataOption
  lazy val jobExecutionMap: JobExecutionMap = standardParams.jobExecutionMap
  lazy val workflowOutputs: CallOutputs = standardParams.workflowOutputs
  override lazy val configurationDescriptor: BackendConfigurationDescriptor = standardParams.configurationDescriptor

  protected val workflowPaths: Option[WorkflowPaths] = initializationDataOption.map {
    _.asInstanceOf[StandardInitializationData].workflowPaths
  }

  protected val ioExecutionContext: MessageDispatcher = context.system.dispatchers.lookup(IoDispatcher)

  override def afterAll(): Future[Unit] = copyCallLogs()

  lazy val logPaths: Seq[Path] = {
    for {
      actualWorkflowPath <- workflowPaths.toSeq
      (backendWorkflowDescriptor, keys) <- jobExecutionMap.toSeq
      key <- keys
      jobPaths = actualWorkflowPath.toJobPaths(key, backendWorkflowDescriptor)
      logPath <- jobPaths.logPaths.values
    } yield logPath
  }

  protected def copyCallLogs(): Future[Unit] = {
    /*
    NOTE: Only using one thread pool slot here to upload all the files for all the calls.
    Using the io-dispatcher defined in application.conf because this might take a while.
    One could also use Future.sequence to flood the dispatcher, or even create a separate jes final call specific thread
    pool for parallel uploads.

    Measure and optimize as necessary. Will likely need retry code at some level as well.
    */
    workflowPaths match {
      case Some(paths) => Future(paths.finalCallLogsPath foreach copyCallLogs)(ioExecutionContext)
      case _ => Future.successful(())
    }
  }

  private def copyCallLogs(callLogsPath: Path): Unit = {
    copyLogs(callLogsPath, logPaths)
  }

  private def copyLogs(callLogsDirPath: Path, logPaths: Seq[Path]): Unit = {
    workflowPaths match {
      case Some(paths) => logPaths.foreach(PathCopier.copy(paths.executionRoot, _, callLogsDirPath))
      case None =>
    }
  }
}
