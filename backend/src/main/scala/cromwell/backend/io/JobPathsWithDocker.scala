package cromwell.backend.io

import java.nio.file.Path

import com.typesafe.config.Config
import cromwell.backend.{BackendWorkflowDescriptor, BackendJobDescriptorKey}
import cromwell.core.path.PathBuilder

class JobPathsWithDocker(val jobKey: BackendJobDescriptorKey,
                         workflowDescriptor: BackendWorkflowDescriptor,
                         config: Config,
                         pathBuilders: List[PathBuilder] = WorkflowPaths.DefaultPathBuilders) extends WorkflowPathsWithDocker(
  workflowDescriptor, config, pathBuilders) with JobPaths {
  import JobPaths._
  
  override lazy val callExecutionRoot = { callRoot.resolve("execution") }
  val callDockerRoot = callPathBuilder(dockerWorkflowRoot, jobKey)
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
        val subpath = p.subpath(executionRoot.getNameCount, p.getNameCount)
        WorkflowPathsWithDocker.DockerRoot.resolve(subpath)
    }
  }
}