package cromwell.backend.google.pipelines.batch

import cromwell.backend.BackendJobDescriptorKey
import cromwell.backend.io.JobPaths

object GcpBatchJobPaths {

}
case class GcpBatchJobPaths(override val workflowPaths: GcpBatchWorkflowPaths, jobKey: BackendJobDescriptorKey, override val isCallCacheCopyAttempt: Boolean = false) extends JobPaths {

  def logBasename = {
    val index = jobKey
      .index
      .map(s => s"-$s")
      .getOrElse("")
    s"${
      jobKey
        .node
        .localName
    }$index"
  }


  override val returnCodeFilename: String = s"$logBasename-rc.txt"
  override def defaultStdoutFilename: String = s"$logBasename-stdout.log"
  override def defaultStderrFilename: String = s"$logBasename-stderr.log"

  override def forCallCacheCopyAttempts: JobPaths = this.copy(isCallCacheCopyAttempt = true)
}
