package cromwell.backend.sfs

import java.nio.file.Path

import cromwell.backend.io.{JobPaths, JobPathsWithDocker, WorkflowPathsBackendInitializationData}
import cromwell.backend.wdl._
import cromwell.backend.{BackendConfigurationDescriptor, BackendInitializationData, BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.CallContext
import cromwell.core.path.PathBuilder
import wdl4s.expression.PureStandardLibraryFunctionsLike
import wdl4s.values.{WdlFile, WdlValue}

import scala.language.postfixOps
import scala.util.{Success, Try}

object SharedFileSystemExpressionFunctions {
  private val LocalFileSystemScheme = "file"

  def isLocalPath(path: Path) = path.toUri.getScheme == SharedFileSystemExpressionFunctions.LocalFileSystemScheme

  def apply(workflowDescriptor: BackendWorkflowDescriptor,
            jobKey: BackendJobDescriptorKey,
            configurationDescriptor: BackendConfigurationDescriptor,
            pathBuilders: List[PathBuilder]): SharedFileSystemExpressionFunctions = {
    val jobPaths = new JobPathsWithDocker(jobKey, workflowDescriptor, configurationDescriptor.backendConfig)
    val callContext = CallContext(
      jobPaths.callExecutionRoot,
      jobPaths.stdout.toString,
      jobPaths.stderr.toString
    )
    new SharedFileSystemExpressionFunctions(pathBuilders, callContext)
  }

  def apply(jobPaths: JobPaths, pathBuilders: List[PathBuilder]): SharedFileSystemExpressionFunctions = {
    val callContext = CallContext(
      jobPaths.callExecutionRoot,
      jobPaths.stdout.toString,
      jobPaths.stderr.toString
    )
    new SharedFileSystemExpressionFunctions(pathBuilders, callContext)
  }

  def apply(workflowDescriptor: BackendWorkflowDescriptor,
            configurationDescriptor: BackendConfigurationDescriptor,
            jobKey: BackendJobDescriptorKey,
            initializationData: Option[BackendInitializationData]) = {
    val jobPaths = new JobPathsWithDocker(jobKey, workflowDescriptor, configurationDescriptor.backendConfig)
    val callContext = CallContext(
      jobPaths.callExecutionRoot,
      jobPaths.stdout.toString,
      jobPaths.stderr.toString
    )

    new SharedFileSystemExpressionFunctions(WorkflowPathsBackendInitializationData.pathBuilders(initializationData), callContext)
  }
}

class SharedFileSystemExpressionFunctions(override val pathBuilders: List[PathBuilder],
                                          context: CallContext
                                 ) extends PureStandardLibraryFunctionsLike with ReadLikeFunctions with WriteFunctions {
  import SharedFileSystemExpressionFunctions._
  import better.files._

  override def writeTempFile(path: String, prefix: String, suffix: String, content: String): String = super[WriteFunctions].writeTempFile(path, prefix, suffix, content)
  override def globPath(glob: String) = context.root.toString
  override def glob(path: String, pattern: String): Seq[String] = {
    File(context.root).glob(s"**/$pattern") map { _.pathAsString } toSeq
  }

  override val writeDirectory = context.root

  override def stdout(params: Seq[Try[WdlValue]]) = Success(WdlFile(context.stdout))
  override def stderr(params: Seq[Try[WdlValue]]) = Success(WdlFile(context.stderr))

  override def postMapping(path: Path) = if (!path.isAbsolute && isLocalPath(path)) context.root.resolve(path) else path
}
