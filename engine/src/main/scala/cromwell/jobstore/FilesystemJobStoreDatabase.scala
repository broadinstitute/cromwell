package cromwell.jobstore

import java.io.File
import java.nio.file.{Files, Paths}

import cromwell.core.ExecutionIndex._
import cromwell.core.WorkflowId
import org.apache.commons.io.FileUtils
import spray.json._

import scala.concurrent.{ExecutionContext, Future}


case object FilesystemJobStoreDatabase extends JobStoreDatabase {

  val separator = File.separator
  val fsroot = Paths.get("JobStore").toAbsolutePath

  implicit class EnhancedJobStoreKey(val jobStoreKey: JobStoreKey) extends AnyVal {
    def fileDir = Paths.get(s"${fsroot.toString}$separator${jobStoreKey.workflowId}").toAbsolutePath
    def filePath = Paths.get(s"$fileDir$separator${jobStoreKey.callFqn}_index${jobStoreKey.index.fromIndex}_attempt${jobStoreKey.attempt}")
  }

  override def writeToDatabase(jobCompletions: Map[JobStoreKey, JobResult], workflowCompletions: List[WorkflowId])(implicit ec: ExecutionContext): Future[Unit] = {
    for {
      _ <- writeJobResultsSerially(jobCompletions)
      _ <- clearWorkflowResultsSerially(workflowCompletions)
    } yield ()
  }

  private def writeJobResultsSerially(jobCompletions: Map[JobStoreKey, JobResult])(implicit ec: ExecutionContext): Future[Unit] = {

    def folderFunction(current: Future[Unit], nextPair: (JobStoreKey, JobResult)): Future[Unit] = {
      val (jobStoreKey, jobResult) = nextPair
      current flatMap { _ => writeJobResult(jobStoreKey, jobResult) }
    }
    jobCompletions.foldLeft(Future.successful(()))(folderFunction)
  }

  private def writeJobResult(jobStoreKey: JobStoreKey, jobResult: JobResult)(implicit ec: ExecutionContext): Future[Unit] = {
    Future {
      import JobResultJsonFormatter._
      val fileDir = jobStoreKey.fileDir
      if (Files.notExists(fileDir)) Files.createDirectories(fileDir)
      FileUtils.writeStringToFile(jobStoreKey.filePath.toFile, jobResult.toJson.toString)
    }
  }

  private def clearWorkflowResultsSerially(workflowCompletions: List[WorkflowId])(implicit ec: ExecutionContext): Future[Unit] = {

    def folderFunction(current: Future[Unit], nextId: WorkflowId): Future[Unit] = current flatMap { _ => clearWorkflowResults(nextId) }

    workflowCompletions.foldLeft(Future.successful(()))(folderFunction)
  }

  private def clearWorkflowResults(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[Unit] = {
    Future {
      val path = Paths.get(s"${fsroot.toString}$separator$workflowId").toAbsolutePath
      if (Files.exists(path)) {
        path.toFile.listFiles foreach { x =>
          Files.delete(x.toPath)
        }
        Files.delete(path)
      }
    }
  }

  override def readJobResult(jobStoreKey: JobStoreKey)(implicit ec: ExecutionContext): Future[Option[JobResult]] = Future {
    import JobResultJsonFormatter._
    val jobFile: File = jobStoreKey.filePath.toFile
    if (jobFile.exists()) Option(FileUtils.readFileToString(jobFile).parseJson.convertTo[JobResult]) else None
  }
}
