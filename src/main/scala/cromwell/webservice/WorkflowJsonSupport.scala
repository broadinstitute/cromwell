package cromwell.webservice

import cromwell.binding.values.WdlValueJsonFormatter._
import cromwell.engine.backend.CallLogs
import spray.json.DefaultJsonProtocol

object WorkflowJsonSupport extends DefaultJsonProtocol {
  implicit val workflowStatusResponseProtocol = jsonFormat2(WorkflowStatusResponse)
  implicit val workflowAbortResponseProtocol = jsonFormat2(WorkflowAbortResponse)
  implicit val workflowSubmitResponseProtocol = jsonFormat2(WorkflowSubmitResponse)
  implicit val workflowOutputResponseProtocol = jsonFormat2(WorkflowOutputResponse)
  implicit val callOutputResponseProtocol = jsonFormat3(CallOutputResponse)
  implicit val callLogsResponseProtocol = jsonFormat2(CallLogs)
  implicit val callStdoutStderrResponse = jsonFormat2(CallStdoutStderrResponse)
}
