package cromwell.backend.io

import java.nio.file.Path

import com.typesafe.config.Config
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}

object JobPaths {
  private val CallPrefix = "call"
  private val ShardPrefix = "shard"
  private val AttemptPrefix = "attempt"
}

class JobPaths(workflowDescriptor: BackendWorkflowDescriptor,
               config: Config,
               jobKey: BackendJobDescriptorKey) extends WorkflowPaths(workflowDescriptor, config) {
  import JobPaths._

  private def callPathBuilder(root: Path) = {
    val callName = jobKey.call.fullyQualifiedName.split('.').last
    val call = s"$CallPrefix-$callName"
    val shard = jobKey.index map { s => s"$ShardPrefix-$s" } getOrElse ""
    val retry = if (jobKey.attempt > 1) s"$AttemptPrefix-${jobKey.attempt}" else ""

    List(call, shard, retry).foldLeft(root)((path, dir) => path.resolve(dir))
  }

  def toDockerPath(path: Path): Path = {
    path.toAbsolutePath match {
      case p if p.startsWith(WorkflowPaths.DockerRoot) => p
      case p =>
        /** For example:
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
        WorkflowPaths.DockerRoot.resolve(subpath)
    }
  }

  val callRoot = callPathBuilder(workflowRoot)
  val callDockerRoot = callPathBuilder(dockerWorkflowRoot)

  val callExecutionRoot = callRoot.resolve("execution")
  val callExecutionDockerRoot = callDockerRoot.resolve("execution")

  val callInputsRoot = callRoot.resolve("inputs")

  val stdout = callExecutionRoot.resolve("stdout")
  val stderr = callExecutionRoot.resolve("stderr")
  val script = callExecutionRoot.resolve("script")
  val returnCode = callExecutionRoot.resolve("rc")
}
