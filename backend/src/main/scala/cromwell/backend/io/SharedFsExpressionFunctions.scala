package cromwell.backend.io

import java.nio.file.{FileSystem, Path}

import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.backend.wdl._
import cromwell.core.CallContext
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.values.{WdlFile, WdlValue}

import scala.language.postfixOps
import scala.util.{Success, Try}

object SharedFsExpressionFunctions {
  private val LocalFSScheme = "file"
  def isLocalPath(path: Path) = path.toUri.getScheme == SharedFsExpressionFunctions.LocalFSScheme
  def apply(workflowDescriptor: BackendWorkflowDescriptor,
            jobKey: BackendJobDescriptorKey,
            configurationDescriptor: BackendConfigurationDescriptor,
            fileSystems: List[FileSystem]): SharedFsExpressionFunctions = {
    val jobPaths = new JobPaths(workflowDescriptor, configurationDescriptor.backendConfig, jobKey)
    val callContext = new CallContext(
      jobPaths.callRoot,
      jobPaths.stdout.toAbsolutePath.toString,
      jobPaths.stderr.toAbsolutePath.toString
    )
    new SharedFsExpressionFunctions(fileSystems, callContext)
  }

  def apply(jobPaths: JobPaths, fileSystems: List[FileSystem]): SharedFsExpressionFunctions = {
    val callContext = new CallContext(
      jobPaths.callRoot,
      jobPaths.stdout.toAbsolutePath.toString,
      jobPaths.stderr.toAbsolutePath.toString
    )
    new SharedFsExpressionFunctions(fileSystems, callContext)
  }

}

class SharedFsExpressionFunctions(override val fileSystems: List[FileSystem],
                                  context: CallContext
                                 ) extends WdlStandardLibraryFunctions with PureFunctions with ReadLikeFunctions with WriteFunctions {
  import SharedFsExpressionFunctions._
  import better.files._

  override def globPath(glob: String) = context.root.toAbsolutePath.toString
  override def glob(path: String, pattern: String): Seq[String] = {
    toPath(path).glob(s"**/$pattern") map { _.path.fullPath } toSeq
  }

  override val writeDirectory = context.root

  override def stdout(params: Seq[Try[WdlValue]]) = Success(WdlFile(context.stdout))
  override def stderr(params: Seq[Try[WdlValue]]) = Success(WdlFile(context.stderr))

  override def postMapping(path: Path) = if (!path.isAbsolute && isLocalPath(path)) context.root.resolve(path) else path
}

