package cromwell.backend.io

import com.typesafe.config.Config
import common.util.StringUtil._
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.path.Path

object JobPathsWithDocker {
  def apply(jobKey: BackendJobDescriptorKey,
            workflowDescriptor: BackendWorkflowDescriptor,
            config: Config) = {
    val workflowPaths = new WorkflowPathsWithDocker(workflowDescriptor, config, WorkflowPaths.DefaultPathBuilders)
    new JobPathsWithDocker(workflowPaths, jobKey)
  }
}

case class JobPathsWithDocker private[io] (override val workflowPaths: WorkflowPathsWithDocker, jobKey: BackendJobDescriptorKey, override val isCallCacheCopyAttempt: Boolean = false) extends JobPaths {
  import JobPaths._

  override lazy val callExecutionRoot = { callRoot.resolve("execution") }
  override def isDocker: Boolean = true
  val callDockerRoot = callPathBuilder(workflowPaths.dockerWorkflowRoot, jobKey, isCallCacheCopyAttempt)
  val callExecutionDockerRoot = callDockerRoot.resolve("execution")
  val callInputsRoot = callRoot.resolve("inputs")
  val callInputsDockerRoot = callDockerRoot.resolve("inputs")

  private lazy val callInputsDockerRootWithSlash = callInputsDockerRoot.pathAsString.ensureSlashed
  private lazy val callExecutionDockerRootWithSlash = callExecutionDockerRoot.pathAsString.ensureSlashed

  override def isInExecution(string: String): Boolean = string.startsWith(callExecutionDockerRootWithSlash)

  override def hostPathFromContainerPath(string: String): Path = {
    callExecutionRoot.resolve(string.stripPrefix(callExecutionDockerRootWithSlash))
  }


  override def hostPathFromContainerInputs(string: String): Path = {
    val stripped = string.stripPrefix(callInputsDockerRootWithSlash)
    callInputsRoot.resolve(stripped)
  }

  def toDockerPath(path: Path): Path = {
    path.toAbsolutePath match {
      case p if p.startsWith(workflowPaths.dockerRoot) => p
      case p =>
        /* For example:
          *
          * p = /abs/path/to/cromwell-executions/three-step/f00ba4/call-ps/stdout.txt
          * localExecutionRoot = /abs/path/to/cromwell-executions
          * subpath = three-step/f00ba4/call-ps/stdout.txt
          *
          * return value = /root/three-step/f00ba4/call-ps/stdout.txt
          *
          * TODO: this assumes that p.startsWith(localExecutionRoot)
          */
        val subpath = p.subpath(workflowPaths.executionRoot.getNameCount, p.getNameCount)
        workflowPaths.dockerRoot.resolve(subpath)
    }
  }

  override def forCallCacheCopyAttempts: JobPaths = this.copy(isCallCacheCopyAttempt = true)
}
