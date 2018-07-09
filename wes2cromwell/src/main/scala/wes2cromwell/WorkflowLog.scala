package wes2cromwell

import spray.json.JsObject

final case class WorkflowLogEntry(
  name: Option[String],
  cmd: Option[Seq[String]],
  start_time: Option[String],
  end_time: Option[String],
  stdout: Option[String],
  stderr: Option[String],
  exit_code: Option[Int]
)

final case class WorkflowLog(
  workflow_id: String,
  request: WorkflowRequest,
  state: WorkflowState,
  workflow_log: Option[WorkflowLogEntry],
  task_logs: Option[Seq[WorkflowLogEntry]],
  outputs: Option[JsObject]
)
