package cromwell.engine.workflow.workflowstore

import java.io.File
import java.nio.file.{Files, Path, Paths}
import java.util.UUID

import cromwell.core.{WorkflowId, WorkflowSourceFiles}
import org.apache.commons.io.FileUtils

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scalaz.NonEmptyList

case object FileSystemWorkflowStore$ extends WorkflowStore {

  val separator = File.separator
  val fsroot= Paths.get("WorkflowStore").toAbsolutePath
  val submittedPath = Paths.get(s"$fsroot${separator}submitted").toAbsolutePath
  val restartablePath = Paths.get(s"$fsroot${separator}restartable").toAbsolutePath
  val runningPath = Paths.get(s"$fsroot${separator}running").toAbsolutePath

  def workflowPath(base: Path, id: WorkflowId): Path = Paths.get(s"$base$separator$id")

  def wdlFilePath(base: Path) = Paths.get(s"${base.toAbsolutePath}${separator}workflow.wdl")
  def inputsFilePath(base: Path) = Paths.get(s"${base.toAbsolutePath}${separator}inputs.json")
  def optionsFilePath(base: Path) = Paths.get(s"${base.toAbsolutePath}${separator}options.json")

  override def initialize(implicit ec: ExecutionContext) = Future {
    if (!Files.exists(submittedPath)) Files.createDirectories(submittedPath)
    if (!Files.exists(restartablePath)) Files.createDirectories(restartablePath)
    if (!Files.exists(runningPath)) Files.createDirectories(runningPath)

    runningPath.toFile.listFiles() foreach { toMove => Files.move(toMove.toPath, Paths.get(s"$restartablePath$separator${toMove.getName}")) }
  }

  /**
    * Adds the requested WorkflowSourceFiles to the store and returns a WorkflowId for each one (in order)
    * for tracking purposes.
    */
  override def add(sources: NonEmptyList[WorkflowSourceFiles])(implicit ec: ExecutionContext): Future[NonEmptyList[WorkflowId]] = {
    val nelOfFutures = sources map { source => Future {
      val id = WorkflowId.randomId()
      val newWorkflowPath = workflowPath(submittedPath, id)
      Files.createDirectory(newWorkflowPath)

      // Write out the files:
      FileUtils.writeStringToFile(wdlFilePath(newWorkflowPath).toFile, source.wdlSource.toString)
      FileUtils.writeStringToFile(inputsFilePath(newWorkflowPath).toFile, source.inputsJson.toString)
      FileUtils.writeStringToFile(optionsFilePath(newWorkflowPath).toFile, source.workflowOptionsJson.toString)
      id
    }}

    // This is horrible because Future.sequence doesn't work with NonEmptyList. Shrug. I assure you it's ok because there's a Nel on the input
    Future.sequence(nelOfFutures.list) map { list => NonEmptyList.nel(list.head, list.tail) }
  }

  override def remove(id: WorkflowId)(implicit ec: ExecutionContext): Future[Boolean] = Future {
    // This is slightly odd but basically: try to delete from all the sub directories. Returns true as soon as any of the
    // delete()s return true.
    List(submittedPath, restartablePath, runningPath) exists { path =>
      val possiblePath = workflowPath(path, id)
      if (Files.exists(possiblePath)) {
        (possiblePath.toFile.listFiles() forall { x: File =>
          x.delete()
        }) && possiblePath.toFile.delete()
      } else {
        false
      }
    }
  }

  override def fetchRunnableWorkflows(n: Int, state: StartableState)(implicit ec: ExecutionContext): Future[List[WorkflowToStart]] = Future {

    val path = if (state == Restartable) restartablePath else submittedPath
    val workflowDirs = path.toFile.listFiles() take n toList

    for {
      wfDir <- workflowDirs
      id = WorkflowId(UUID.fromString(wfDir.getName))
      wdlFile = wdlFilePath(wfDir.toPath).toFile
      wdl = FileUtils.readFileToString(wdlFile)
      inputs = FileUtils.readFileToString(inputsFilePath(wfDir.toPath).toFile)
      options = FileUtils.readFileToString(optionsFilePath(wfDir.toPath).toFile)
      _ = Files.move(wfDir.toPath, Paths.get(s"$runningPath$separator${wfDir.getName}"))
    } yield WorkflowToStart(id, WorkflowSourceFiles(wdl, inputs, options), state)
  }
}
