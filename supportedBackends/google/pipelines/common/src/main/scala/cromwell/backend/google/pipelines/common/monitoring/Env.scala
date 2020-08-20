package cromwell.backend.google.pipelines.common.monitoring

import cromwell.backend.BackendJobDescriptor

object Env {
  /**
    * Name of an environmental variable
    */
  val WorkflowId = "WORKFLOW_ID"
  val TaskCallName = "TASK_CALL_NAME"
  val TaskCallIndex = "TASK_CALL_INDEX"
  val TaskCallAttempt = "TASK_CALL_ATTEMPT"
  val DiskMounts = "DISK_MOUNTS"

  def monitoringImageEnvironment(jobDescriptor: BackendJobDescriptor)
                                (mountPaths: List[String]): Map[String, String] =
    Map(
      Env.WorkflowId -> jobDescriptor.workflowDescriptor.id.toString,
      Env.TaskCallName -> jobDescriptor.taskCall.localName,
      Env.TaskCallIndex -> jobDescriptor.key.index.map(_.toString).getOrElse("NA"),
      Env.TaskCallAttempt -> jobDescriptor.key.attempt.toString,
      Env.DiskMounts -> mountPaths.mkString(" "),
    )
}
