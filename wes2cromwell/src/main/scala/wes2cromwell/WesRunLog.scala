package wes2cromwell

import akka.http.scaladsl.marshallers.sprayjson.SprayJsonSupport
import spray.json.{DefaultJsonProtocol, JsObject, JsonFormat}

final case class WesLog(
  name: Option[String],
  cmd: Option[Seq[String]],
  start_time: Option[String],
  end_time: Option[String],
  stdout: Option[String],
  stderr: Option[String],
  exit_code: Option[Int]
)

final case class WesRunLog(run_id: String,
                           state: WesRunState,
                           run_log: Option[WesLog],
                           task_logs: Option[Seq[WesLog]],
                           outputs: Option[JsObject])

object WorkflowLogJsonSupport extends SprayJsonSupport with DefaultJsonProtocol {
  implicit val workflowLogEntryFormat: JsonFormat[WesLog] = jsonFormat7(WesLog)
  implicit val workflowLogFormat: JsonFormat[WesRunLog] = jsonFormat5(WesRunLog)
}