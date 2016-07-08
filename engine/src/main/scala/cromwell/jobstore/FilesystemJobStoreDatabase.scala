package cromwell.jobstore

import java.io.File
import java.nio.file.{Files, Path, Paths}

import cromwell.core.ExecutionIndex._
import cromwell.core.{JobOutputs, WorkflowId}
import spray.json._
import wdl4s.WdlExpression
import wdl4s.types.WdlArrayType
import wdl4s.values._
import cromwell.webservice.WdlValueJsonFormatter._
import org.apache.commons.io.FileUtils

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try


case object FilesystemJobStoreDatabase extends JobStoreDatabase {

  val separator = File.separator
  val fsroot = Paths.get("JobStore").toAbsolutePath

  override def writeToDatabase(jobCompletions: Map[JobStoreKey, JobResult], workflowCompletions: List[WorkflowId])(implicit ec: ExecutionContext): Future[Unit] = {
    for {
      _ <- writeJobResultsSerially(jobCompletions)
      _ <- writeWorkflowResultsSerially(workflowCompletions)
    } yield ()
  }

  private def writeJobResultsSerially(jobCompletions: Map[JobStoreKey, JobResult])(implicit ec: ExecutionContext): Future[Unit] = {

    def folderFunction(current: Future[Unit], nextPair: (JobStoreKey, JobResult)): Future[Unit] = {
      val (jobStoreKey, jobResult) = nextPair
      current map { _ => writeJobResult(jobStoreKey, jobResult) }
    }
    jobCompletions.foldLeft(Future.successful(()))(folderFunction)
  }

  private def writeJobResult(jobStoreKey: JobStoreKey, jobResult: JobResult)(implicit ec: ExecutionContext): Future[Unit] = {
    Future {
      import JobResultJsonFormatter._
      val fileDir = Paths.get(s"${fsroot.toString}$separator${jobStoreKey.workflowId}").toAbsolutePath
      if (Files.notExists(fileDir)) Files.createDirectories(fileDir)
      val filePath = Paths.get(s"$fileDir$separator${jobStoreKey.callFqn}_index${jobStoreKey.index.fromIndex}_attempt${jobStoreKey.attempt}")

      FileUtils.writeStringToFile(filePath.toFile, jobResult.toJson.toString)
    }
  }

  private def writeWorkflowResultsSerially(workflowCompletions: List[WorkflowId])(implicit ec: ExecutionContext): Future[Unit] = {

    def folderFunction(current: Future[Unit], nextId: WorkflowId): Future[Unit] = current map { _ => clearWorkflowResults(nextId) }

    workflowCompletions.foldLeft(Future.successful(()))(folderFunction)
  }

  private def clearWorkflowResults(workflowId: WorkflowId)(implicit ec: ExecutionContext): Future[Unit] = {
    Future {
      val path = Paths.get(s"${fsroot.toString}$separator$workflowId").toAbsolutePath
      if (Files.exists(path)) {
        //Files.delete(path)
      }
    }
  }
}
