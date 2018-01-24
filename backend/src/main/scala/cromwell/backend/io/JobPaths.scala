package cromwell.backend.io

import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.path.Path
import cromwell.core.{CallContext, JobKey}
import cromwell.services.metadata.CallMetadataKeys
import wom.expression.WomExpression

object JobPaths {
  val CallPrefix = "call"
  val ShardPrefix = "shard"
  val AttemptPrefix = "attempt"
  val ScriptPathKey = "script"
  val StdoutPathKey = "stdout"
  val StdErrPathKey = "stderr"
  val ReturnCodePathKey = "returnCode"
  val CallRootPathKey = "callRootPath"
  val DockerCidPathKey = "dockerCidPath"

  def callPathBuilder(root: Path, jobKey: JobKey) = {
    val callName = jobKey.node.localName
    val call = s"$CallPrefix-$callName"
    val shard = jobKey.index map { s => s"$ShardPrefix-$s" } getOrElse ""
    val retry = if (jobKey.attempt > 1) s"$AttemptPrefix-${jobKey.attempt}" else ""

    List(call, shard, retry).foldLeft(root)((path, dir) => path.resolve(dir))
  }
}

trait JobPaths {
  import JobPaths._

  def workflowPaths: WorkflowPaths
  def returnCodeFilename: String = "rc"
  def defaultStdoutFilename = "stdout"
  def defaultStderrFilename = "stderr"

  final def stdinFilenameExpression: Option[WomExpression] = jobKey.node.callable.stdinRedirection
  final def stdoutFilename: String = jobKey.node.callable.stdoutRedirection.getOrElse(defaultStdoutFilename)
  final def stderrFilename: String = jobKey.node.callable.stderrRedirection.getOrElse(defaultStderrFilename)
  def scriptFilename: String = "script"
  def dockerCidFilename: String = "docker_cid"

  def jobKey: BackendJobDescriptorKey
  lazy val callRoot = callPathBuilder(workflowPaths.workflowRoot, jobKey)
  lazy val callExecutionRoot = callRoot
  lazy val stdout = callExecutionRoot.resolve(stdoutFilename)
  lazy val stderr = callExecutionRoot.resolve(stderrFilename)
  lazy val script = callExecutionRoot.resolve(scriptFilename)
  lazy val dockerCid = callExecutionRoot.resolve(dockerCidFilename)
  lazy val returnCode = callExecutionRoot.resolve(returnCodeFilename)

  private lazy val commonMetadataPaths: Map[String, Path] = Map(
    CallMetadataKeys.CallRoot -> callRoot,
    CallMetadataKeys.Stdout -> stdout,
    CallMetadataKeys.Stderr -> stderr
  )

  private lazy val commonDetritusPaths: Map[String, Path] = Map(
    JobPaths.CallRootPathKey -> callRoot,
    JobPaths.ScriptPathKey -> script,
    JobPaths.StdoutPathKey -> stdout,
    JobPaths.StdErrPathKey -> stderr,
    JobPaths.ReturnCodePathKey -> returnCode
  )

  private lazy val commonLogPaths: Map[String, Path] = Map(
    JobPaths.StdoutPathKey -> stdout,
    JobPaths.StdErrPathKey -> stderr
  )

  protected lazy val customMetadataPaths: Map[String, Path] = Map.empty
  protected lazy val customDetritusPaths: Map[String, Path] = Map.empty
  protected lazy val customLogPaths: Map[String, Path] = Map.empty

  lazy val metadataPaths = commonMetadataPaths ++ customMetadataPaths
  lazy val detritusPaths = commonDetritusPaths ++ customDetritusPaths
  lazy val logPaths = commonLogPaths ++ customLogPaths

  lazy val callContext = CallContext(callExecutionRoot, stdout.pathAsString, stderr.pathAsString)
}
