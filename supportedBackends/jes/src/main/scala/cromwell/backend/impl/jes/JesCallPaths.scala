package cromwell.backend.impl.jes

import java.nio.file.Path

import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.CallContext

object JesCallPaths {
  def apply(jobKey: BackendJobDescriptorKey, workflowDescriptor: BackendWorkflowDescriptor, jesConfiguration: JesConfiguration): JesCallPaths = {
    new JesCallPaths(jobKey, workflowDescriptor, jesConfiguration)
  }
}

class JesCallPaths(jobKey: BackendJobDescriptorKey, workflowDescriptor: BackendWorkflowDescriptor, jesConfiguration: JesConfiguration
                  ) extends JesWorkflowPaths(workflowDescriptor, jesConfiguration) {

  val CallPrefix = "call"
  val ShardPrefix = "shard"
  val AttemptPrefix = "attempt"

  val jesLogBasename = {
    val index = jobKey.index.map(s => s"-$s").getOrElse("")
    s"${jobKey.scope.unqualifiedName}$index"
  }

  val callRootPath: Path = {
    val callName = jobKey.call.fullyQualifiedName.split('.').last
    val call = s"$CallPrefix-$callName"
    val shard = jobKey.index map { s => s"$ShardPrefix-$s" } getOrElse ""
    val retry = if (jobKey.attempt > 1) s"$AttemptPrefix-${jobKey.attempt}" else ""

    List(call, shard, retry).foldLeft(workflowRootPath)((path, dir) => path.resolve(dir))
  }

  val returnCodeFilename: String = s"$jesLogBasename-rc.txt"
  val stdoutFilename: String = s"$jesLogBasename-stdout.log"
  val stderrFilename: String = s"$jesLogBasename-stderr.log"
  val jesLogFilename: String = s"$jesLogBasename.log"

  lazy val returnCodePath: Path = callRootPath.resolve(returnCodeFilename)
  lazy val stdoutPath: Path = callRootPath.resolve(stdoutFilename)
  lazy val stderrPath: Path = callRootPath.resolve(stderrFilename)
  lazy val jesLogPath: Path = callRootPath.resolve(jesLogFilename)
  lazy val callContext = new CallContext(callRootPath, stdoutFilename, stderrFilename)
}
