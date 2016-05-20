package cromwell.backend.impl.jes

import java.nio.file.Path

import com.typesafe.config.Config
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.CallContext
import cromwell.filesystems.gcs.GcsFileSystem

object JesCallPaths {
  def apply(jobKey: BackendJobDescriptorKey, gcsFileSystem: GcsFileSystem, workflowDescriptor: BackendWorkflowDescriptor, backendConfig: Config): JesCallPaths = {
    val GcsRootOptionKey = "jes_gcs_root"
    val CallPrefix = "call"
    val ShardPrefix = "shard"
    val AttemptPrefix = "attempt"

    def jesLogBasename = {
      val index = jobKey.index.map(s => s"-$s").getOrElse("")
      s"${jobKey.scope.unqualifiedName}$index"
    }

    val rootPath: Path =
      gcsFileSystem.getPath(workflowDescriptor.workflowOptions.getOrElse(GcsRootOptionKey, backendConfig.getString("root")))

    val workflowRootPath: Path = rootPath.resolve(workflowDescriptor.workflowNamespace.workflow.unqualifiedName)
      .resolve(workflowDescriptor.id.toString)

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
    JesCallPaths(rootPath, workflowRootPath, callRootPath, stdoutFilename, stderrFilename, jesLogFilename, returnCodeFilename)
  }
}

case class JesCallPaths private(rootPath: Path, workflowRootPath: Path, callRootPath: Path, stdoutFilename: String, stderrFilename: String, jesLogFilename: String, returnCodeFilename: String) {
  lazy val returnCodePath: Path = callRootPath.resolve(returnCodeFilename)
  lazy val stdoutPath: Path = callRootPath.resolve(stdoutFilename)
  lazy val stderrPath: Path = callRootPath.resolve(stderrFilename)
  lazy val callContext = new CallContext(callRootPath, stdoutFilename, stderrFilename)
}