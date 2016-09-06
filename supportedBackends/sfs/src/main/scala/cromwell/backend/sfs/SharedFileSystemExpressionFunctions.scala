package cromwell.backend.sfs

import java.nio.file.{FileSystem, Path}

import cromwell.backend.io.{JobPaths, WorkflowPathsBackendInitializationData}
import cromwell.backend.wdl._
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.CallContext
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.values.{WdlFile, WdlValue}

import scala.language.postfixOps
import scala.util.{Success, Try}

object SharedFileSystemExpressionFunctions {
  private val LocalFileSystemScheme = "file"

  def isLocalPath(path: Path) = path.toUri.getScheme == SharedFileSystemExpressionFunctions.LocalFileSystemScheme

  def apply(workflowDescriptor: BackendWorkflowDescriptor,
            jobKey: BackendJobDescriptorKey,
            configurationDescriptor: BackendConfigurationDescriptor,
            fileSystems: List[FileSystem]): SharedFileSystemExpressionFunctions = {
    val jobPaths = new JobPaths(workflowDescriptor, configurationDescriptor.backendConfig, jobKey)
    val callContext = CallContext(
      jobPaths.callRoot,
      jobPaths.stdout.toString,
      jobPaths.stderr.toString
    )
    new SharedFileSystemExpressionFunctions(fileSystems, callContext)
  }

  def apply(jobPaths: JobPaths, fileSystems: List[FileSystem]): SharedFileSystemExpressionFunctions = {
    val callContext = CallContext(
      jobPaths.callRoot,
      jobPaths.stdout.toString,
      jobPaths.stderr.toString
    )
    new SharedFileSystemExpressionFunctions(fileSystems, callContext)
  }

  def apply(workflowDescriptor: BackendWorkflowDescriptor,
            configurationDescriptor: BackendConfigurationDescriptor,
            jobKey: BackendJobDescriptorKey,
            initializationData: Option[BackendInitializationData]) = {
    val jobPaths = new JobPaths(workflowDescriptor, configurationDescriptor.backendConfig, jobKey)
    val callContext = CallContext(
      jobPaths.callRoot,
      jobPaths.stdout.toString,
      jobPaths.stderr.toString
    )

    new SharedFileSystemExpressionFunctions(WorkflowPathsBackendInitializationData.fileSystems(initializationData), callContext)
  }
}

class SharedFileSystemExpressionFunctions(override val fileSystems: List[FileSystem],
                                          context: CallContext
                                 ) extends WdlStandardLibraryFunctions with PureFunctions with ReadLikeFunctions with WriteFunctions {
  import SharedFileSystemExpressionFunctions._
  import better.files._

  override def globPath(glob: String) = context.root.toString
  override def glob(path: String, pattern: String): Seq[String] = {
    File(toPath(path)).glob(s"**/$pattern") map { _.pathAsString } toSeq
  }

  override val writeDirectory = context.root

  override def stdout(params: Seq[Try[WdlValue]]) = Success(WdlFile(context.stdout))
  override def stderr(params: Seq[Try[WdlValue]]) = Success(WdlFile(context.stderr))

  override def postMapping(path: Path) = if (!path.isAbsolute && isLocalPath(path)) context.root.resolve(path) else path
}
