package cromwell.backend.io

import common.util.StringUtil._
import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.path.Path
import cromwell.core.{CallContext, JobKey, StandardPaths}
import cromwell.services.metadata.CallMetadataKeys

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
  def doubleMemoryRetryRCFilename: String = "double_memory_retry_rc"
  def defaultStdoutFilename = "stdout"
  def defaultStderrFilename = "stderr"
  def isDocker: Boolean = false

  // In this non-Docker version of `JobPaths` there is no distinction between host and container roots so this is
  // just called 'rootWithSlash'.
  private lazy val rootWithSlash = callExecutionRoot.pathAsString.ensureSlashed

  def isInExecution(string: String): Boolean = string.startsWith(rootWithSlash)

  /**
    * Return a host path corresponding to the specified container path.
    */
  def hostPathFromContainerPath(string: String): Path = {
    // No container here, just return a Path of the absolute path to the file.
    callExecutionRoot.resolve(string.stripPrefix(rootWithSlash))
  }

  def hostPathFromContainerInputs(string: String): Path =
    // No container here, just return a Path of the absolute path to the file.
    callExecutionRoot.resolve(string.stripPrefix(rootWithSlash))


  def scriptFilename: String = "script"
  def dockerCidFilename: String = "docker_cid"

  def jobKey: BackendJobDescriptorKey
  lazy val callRoot = callPathBuilder(workflowPaths.workflowRoot, jobKey)
  lazy val callExecutionRoot = callRoot

  // Use default stdout and stderr names by default. This StandardPaths `var` may be reassigned later to
  // enable dynamic standard output and error file names for languages like CWL that support this feature.
  var standardPaths: StandardPaths = StandardPaths(
    output = callExecutionRoot.resolve(defaultStdoutFilename),
    error = callExecutionRoot.resolve(defaultStderrFilename)
  )

  lazy val script = callExecutionRoot.resolve(scriptFilename)
  lazy val dockerCid = callExecutionRoot.resolve(dockerCidFilename)
  lazy val returnCode = callExecutionRoot.resolve(returnCodeFilename)
  lazy val doubleMemoryRetryRC = callExecutionRoot.resolve(doubleMemoryRetryRCFilename)

  // This is a `def` because `standardPaths` is a `var` that may be reassigned during the calculation of
  // standard output and error file names.
  def standardOutputAndErrorPaths: Map[String, Path] = Map(
    CallMetadataKeys.Stdout -> standardPaths.output,
    CallMetadataKeys.Stderr -> standardPaths.error
  )

  private lazy val commonMetadataPaths: Map[String, Path] =
    standardOutputAndErrorPaths + (CallMetadataKeys.CallRoot -> callRoot)

  // This is a `def` because `standardPaths` is a `var` that may be reassigned during the calculation of
  // standard output and error file names.
  private def commonDetritusPaths: Map[String, Path] = Map(
    JobPaths.CallRootPathKey -> callRoot,
    JobPaths.ScriptPathKey -> script,
    JobPaths.StdoutPathKey -> standardPaths.output,
    JobPaths.StdErrPathKey -> standardPaths.error,
    JobPaths.ReturnCodePathKey -> returnCode
  )

  private lazy val commonLogPaths: Map[String, Path] = Map(
    JobPaths.StdoutPathKey -> standardPaths.output,
    JobPaths.StdErrPathKey -> standardPaths.error
  )

  protected lazy val customMetadataPaths: Map[String, Path] = Map.empty
  protected lazy val customDetritusPaths: Map[String, Path] = Map.empty
  protected lazy val customLogPaths: Map[String, Path] = Map.empty

  lazy val metadataPaths = commonMetadataPaths ++ customMetadataPaths
  def detritusPaths = commonDetritusPaths ++ customDetritusPaths
  lazy val logPaths = commonLogPaths ++ customLogPaths

  lazy val callContext = CallContext(callExecutionRoot, standardPaths, isDocker)
}
