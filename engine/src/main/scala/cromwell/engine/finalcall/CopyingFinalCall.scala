package cromwell.engine.finalcall

import java.nio.file.Path

import cromwell.engine.backend.WorkflowDescriptor
import cromwell.webservice.WorkflowMetadataResponse
import wdl4s._
import wdl4s.values.{WdlFile, WdlSingleFile, WdlValue}

trait CopyingFinalCall extends FinalCall {
  def createFinalCallCopies(workflow: WorkflowDescriptor, metadata: WorkflowMetadataResponse): Seq[FinalCallCopy]

  lazy val wdlSource = FinalCallCopyHack.wdlSource(finalCallName) // TODO: Don't hack

  final override def createCallInputs(workflow: WorkflowDescriptor, metadata: WorkflowMetadataResponse): CallInputs = {
    val copies = createFinalCallCopies(workflow, metadata)
    FinalCallCopyHack.toCallInputs(copies) // TODO: Don't hack
  }
}

case class FinalCallCopy(source: Path, destination: Path)

object CopyingFinalCall {
  def copyToDir(workflow: WorkflowDescriptor, destDirectory: String)(wdlValue: WdlValue): Seq[FinalCallCopy] = {
    // All outputs should have wdl values at this point, if they don't there's nothing we can do here
    for {
      wdlFile <- wdlValue collectAsSeq { case f: WdlSingleFile => f }
    } yield copyWdlFile(workflow, wdlFile, destDirectory)
  }

  private def copyWdlFile(workflow: WorkflowDescriptor,
                          file: WdlFile,
                          destDirectory: String): FinalCallCopy = {
    val srcPath = workflow.toFileSystemPath(file.valueString)
    val workflowRootPath = workflow.workflowRootPath.toAbsolutePath
    val srcSubPath = srcPath.subpath(workflowRootPath.getNameCount, srcPath.getNameCount)
    val destRootPath = workflow.toFileSystemPath(destDirectory).toAbsolutePath
    val destWorkflowPath = destRootPath.resolve(workflow.relativeWorkflowRootPath)
    // UnixPath.resolve(NioGcsPath) seems to be returning a null pointer. TODO: Add a test to confirm
    val destPath = destWorkflowPath.resolve(srcSubPath.toString)
    FinalCallCopy(srcPath, destPath)
  }
}
