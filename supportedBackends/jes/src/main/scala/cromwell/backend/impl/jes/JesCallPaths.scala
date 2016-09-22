package cromwell.backend.impl.jes

import java.nio.file.Path

import cromwell.backend.impl.jes.authentication.JesCredentials
import cromwell.backend.io.JobPaths
import cromwell.backend.io.JobPaths._
import cromwell.backend.{BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.CallContext
import cromwell.services.metadata.CallMetadataKeys

import scala.concurrent.ExecutionContext

object JesCallPaths {
  def apply(jobKey: BackendJobDescriptorKey, workflowDescriptor: BackendWorkflowDescriptor,
            jesConfiguration: JesConfiguration,
            credentials: JesCredentials)(implicit ec: ExecutionContext): JesCallPaths = {
    new JesCallPaths(jobKey, workflowDescriptor, jesConfiguration, credentials)
  }

  val JesLogPathKey = "jesLog"
  val GcsExecPathKey = "gcsExec"
}

class JesCallPaths(jobKey: BackendJobDescriptorKey, workflowDescriptor: BackendWorkflowDescriptor,
                   jesConfiguration: JesConfiguration,
                   credentials: JesCredentials)(implicit ec: ExecutionContext) extends
  JesWorkflowPaths(workflowDescriptor, jesConfiguration, credentials)(ec) {

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
  val gcsExecFilename: String = "exec.sh"

  lazy val returnCodePath: Path = callRootPath.resolve(returnCodeFilename)
  lazy val stdoutPath: Path = callRootPath.resolve(stdoutFilename)
  lazy val stderrPath: Path = callRootPath.resolve(stderrFilename)
  lazy val jesLogPath: Path = callRootPath.resolve(jesLogFilename)
  lazy val gcsExecPath: Path = callRootPath.resolve(gcsExecFilename)
  lazy val callContext = CallContext(callRootPath, stdoutFilename, stderrFilename)

  /*
  TODO: Move various monitoring files path generation here.

  "/cromwell_root" is a well known path, called in the regular JobPaths callDockerRoot.
  This JesCallPaths should know about that root, and be able to create the monitoring file paths.
  Instead of the AsyncActor creating the paths, the paths could then be shared with the CachingActor.

  Those monitoring paths could then be returned by metadataFiles and detritusFiles.
   */

  lazy val metadataPaths: Map[String, Path] = Map(
    CallMetadataKeys.CallRoot -> callRootPath,
    CallMetadataKeys.Stdout -> stdoutPath,
    CallMetadataKeys.Stderr -> stderrPath,
    CallMetadataKeys.BackendLogsPrefix + ":log" -> jesLogPath
  )

  lazy val detritusPaths: Map[String, Path] = Map(
    JobPaths.CallRootPathKey -> callRootPath,
    JesCallPaths.GcsExecPathKey -> gcsExecPath,
    JesCallPaths.JesLogPathKey -> jesLogPath,
    JobPaths.StdoutPathKey -> stdoutPath,
    JobPaths.StdErrPathKey -> stderrPath,
    JobPaths.ReturnCodePathKey -> returnCodePath
  )
}
