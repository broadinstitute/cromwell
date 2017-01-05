package cromwell.backend.sfs

import java.nio.file.Path

import cromwell.backend.io._
import cromwell.backend.wdl._
import cromwell.backend._
import cromwell.core.CallContext
import cromwell.core.path.PathBuilder
import wdl4s.expression.PureStandardLibraryFunctionsLike
import wdl4s.values._

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
                                 ) extends PureStandardLibraryFunctionsLike with ReadLikeFunctions with WriteFunctions with GlobFunctions {
  import SharedFileSystemExpressionFunctions._

  def callContext: CallContext = context

  override def writeTempFile(path: String, prefix: String, suffix: String, content: String): String = super[WriteFunctions].writeTempFile(path, prefix, suffix, content)

  override val writeDirectory = context.root

  override def stdout(params: Seq[Try[WdlValue]]) = Success(WdlFile(context.stdout))
  override def stderr(params: Seq[Try[WdlValue]]) = Success(WdlFile(context.stderr))

  override def postMapping(path: Path) = if (!path.isAbsolute && isLocalPath(path)) context.root.resolve(path) else path
}
