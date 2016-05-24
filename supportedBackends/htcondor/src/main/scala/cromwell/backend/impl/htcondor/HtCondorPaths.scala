package cromwell.backend.impl.htcondor


import java.nio.file.{Path, Paths}

import com.typesafe.config.Config
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}

object WorkflowPaths{
  val DockerRoot = Paths.get("/root")
}

class WorkflowPaths(workflowDescriptor: BackendWorkflowDescriptor, config: Config) {
  // TODO will need update when configuration is rearranged
  protected val executionRoot = Paths.get(config.getString("root"))

  private def workflowPathBuilder(root: Path) = {
    root.resolve(workflowDescriptor.workflowNamespace.workflow.unqualifiedName)
      .resolve(workflowDescriptor.id.toString)
  }

  lazy val workflowRoot = workflowPathBuilder(executionRoot)
  lazy val dockerWorkflowRoot = workflowPathBuilder(WorkflowPaths.DockerRoot)
}

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
        val subpath = p.subpath(executionRoot.toAbsolutePath.getNameCount, p.getNameCount)
        WorkflowPaths.DockerRoot.resolve(subpath)
    }
  }

  private val jobName = jobKey.call.unqualifiedName

  val callRoot = callPathBuilder(workflowRoot)
  val callDockerRoot = callPathBuilder(dockerWorkflowRoot)

  val stdout = callRoot.resolve("stdout")
  val stderr = callRoot.resolve("stderr")
  val script = callRoot.resolve(s"$jobName-script.sh")
  val returnCode = callRoot.resolve("rc")
  val submitFile = callRoot.resolve("submitfile")
  val submitFileStderr = callRoot.resolve("submitfile.stderr")
  val submitFileStdout = callRoot.resolve("submitfile.stdout")
  val htcondorLog = callRoot.resolve(s"$jobName.log")
}
