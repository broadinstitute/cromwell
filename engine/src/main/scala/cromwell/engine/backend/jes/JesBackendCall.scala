package cromwell.engine.backend.jes

import java.nio.file.Paths

import com.google.api.client.googleapis.json.GoogleJsonResponseException
import com.typesafe.scalalogging.LazyLogging
import cromwell.engine.backend.jes.Run.TerminalRunStatus
import cromwell.engine.backend.jes.authentication.ProductionJesAuthentication
import cromwell.engine.backend.{BackendCall, CallLogs, JobKey, _}
import cromwell.engine.io.gcs.GcsPath
import cromwell.engine.AbortRegistrationFunction
import wdl4s.values._

import scala.concurrent.{ExecutionContext, Future}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}


class JesBackendCall(val backend: JesBackend,
                     val jobDescriptor: BackendCallJobDescriptor,
                     val callAbortRegistrationFunction: Option[AbortRegistrationFunction])
  extends BackendCall with ProductionJesAuthentication with LazyLogging {

  import JesBackend._

  // TODO: Assuming that runtimeAttributes.disks always has a 'local-disk'
  lazy val workingDisk = runtimeAttributes.disks.find(_.name == JesWorkingDisk.Name).get

  def jesCommandLine = s"/bin/bash ${cmdInput.containerPath.toAbsolutePath.toString}"

  lazy val callGcsPath = key.callRootPathWithBaseRoot(workflowDescriptor, backend.rootPath(workflowDescriptor.workflowOptions))
  //TODO: Switch to NioGcsPath
  lazy val callDir = GcsPath(callGcsPath.toString)
  lazy val gcsExecPath = GcsPath(callGcsPath + "/" + JesExecScript)
  lazy val defaultMonitoringOutputPath = callGcsPath + "/" + JesMonitoringLogFile

  lazy val returnCodeFilename = jesReturnCodeFilename(key)
  lazy val jesStdoutGcsPath = s"$callGcsPath/${jesLogStdoutFilename(key)}"
  lazy val jesStderrGcsPath = s"$callGcsPath/${jesLogStderrFilename(key)}"
  lazy val jesLogGcsPath = s"$callGcsPath/${jesLogFilename(key)}"
  lazy val returnCodeGcsPath = s"$callGcsPath/$returnCodeFilename"

  lazy val rcJesOutput = JesFileOutput(returnCodeFilename, returnCodeGcsPath, Paths.get(returnCodeFilename), workingDisk)
  lazy val cmdInput = JesFileInput(ExecParamName, gcsExecPath.toString, Paths.get(JesExecScript), workingDisk)

  def standardParameters = Seq(rcJesOutput)

  def downloadRcFile = authenticateAsUser(workflowDescriptor) { storage => Try(storage.readFile(returnCodeGcsPath)) }

  /**
    * Determines the maximum number of times a call can be started with a Preemptible VM.
    * TODO: Use configuration as a way to set this globally.
    * Currently workflow options act as default for runtime attributes, configuration could do the same for workflow options.
    */
  lazy val maxPreemption = runtimeAttributes.preemptible

  lazy val preemptible = key.attempt <= maxPreemption

  /**
   * Determine the output directory for the files matching a particular glob.
   */
  def globOutputPath(glob: String) = s"$callGcsPath/glob-${glob.md5Sum}/"

  override def execute(implicit ec: ExecutionContext) = backend.execute(this)

  override def resume(jobKey: JobKey)(implicit ec: ExecutionContext) = backend.resume(this, jobKey)

  override def useCachedCall(avoidedTo: BackendCall)(implicit ec: ExecutionContext): Future[ExecutionHandle] =
    backend.useCachedCall(avoidedTo.asInstanceOf[JesBackendCall], this)

  override def stdoutStderr: CallLogs = {
    CallLogs(
      stdout = WdlFile(jesStdoutGcsPath),
      stderr = WdlFile(jesStderrGcsPath),
      Option(Map("log" -> WdlFile(jesLogGcsPath)))
    )
  }
}
