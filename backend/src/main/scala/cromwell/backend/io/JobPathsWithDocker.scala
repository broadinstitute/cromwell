package cromwell.backend.io

import com.typesafe.config.Config
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.path.{Path, PathBuilder}

object JobPathsWithDocker {
  def apply(jobKey: BackendJobDescriptorKey,
  workflowDescriptor: BackendWorkflowDescriptor,
  config: Config,
  pathBuilders: List[PathBuilder] = WorkflowPaths.DefaultPathBuilders) = {
    val workflowPaths = new WorkflowPathsWithDocker(workflowDescriptor, config, pathBuilders)
    new JobPathsWithDocker(workflowPaths, jobKey)
  }
}
class JobPathsWithDocker(override val workflowPaths: WorkflowPathsWithDocker, val jobKey: BackendJobDescriptorKey) extends JobPaths {
  import JobPaths._
  
  override lazy val callExecutionRoot = { callRoot.resolve("execution") }
  val callDockerRoot = callPathBuilder(workflowPaths.dockerWorkflowRoot, jobKey)
  val callExecutionDockerRoot = callDockerRoot.resolve("execution")
  val callInputsRoot = callRoot.resolve("inputs")

  def toDockerPath(path: Path): Path = {
    path.toAbsolutePath match {
      case p if p.startsWith(WorkflowPathsWithDocker.DockerRoot) => p
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
        WorkflowPathsWithDocker.DockerRoot.resolve(subpath)
    }
  }
}
