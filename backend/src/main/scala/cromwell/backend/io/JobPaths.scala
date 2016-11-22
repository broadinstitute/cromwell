package cromwell.backend.io

import java.nio.file.Path

import cromwell.core.JobKey
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

  def callPathBuilder(root: Path, jobKey: JobKey) = {
    val callName = jobKey.scope.unqualifiedName
    val call = s"$CallPrefix-$callName"
    val shard = jobKey.index map { s => s"$ShardPrefix-$s" } getOrElse ""
    val retry = if (jobKey.attempt > 1) s"$AttemptPrefix-${jobKey.attempt}" else ""

    List(call, shard, retry).foldLeft(root)((path, dir) => path.resolve(dir))
  }
}

trait JobPaths { this: WorkflowPaths =>
  import JobPaths._

  def returnCodeFilename: String = "rc"
  def stdoutFilename: String = "stdout"
  def stderrFilename: String = "stderr"
  def scriptFilename: String = "script"
  
  def jobKey: JobKey
  lazy val callRoot = callPathBuilder(workflowRoot, jobKey)
  lazy val callExecutionRoot = callRoot
  lazy val stdout = callExecutionRoot.resolve(stdoutFilename)
  lazy val stderr = callExecutionRoot.resolve(stderrFilename)
  lazy val script = callExecutionRoot.resolve(scriptFilename)
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
  
  protected lazy val customMetadataPaths: Map[String, Path] = Map.empty
  protected lazy val customDetritusPaths: Map[String, Path] = Map.empty
  
  lazy val metadataPaths = commonMetadataPaths ++ customMetadataPaths
  lazy val detritusPaths = commonDetritusPaths ++ customDetritusPaths
}
