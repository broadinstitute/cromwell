package cromwell.webservice.routes.wes

import cromwell.webservice.routes.wes.WesState.WesState
import spray.json.JsObject


final case class WesLog(name: Option[String],
                        cmd: Option[Seq[String]],
                        start_time: Option[String],
                        end_time: Option[String],
                        stdout: Option[String],
                        stderr: Option[String],
                        exit_code: Option[Int]
                       )

final case class WesRunRequest(workflow_params: Option[JsObject],
                               workflow_type: String,
                               workflow_type_version: String,
                               tags: Option[JsObject],
                               workflow_engine_parameters: Option[JsObject],
                               workflow_url: Option[String]
                              )

final case class WesRunLog(run_id: String,
                           request: WesRunRequest,
                           state: WesState,
                           run_log: Option[WesLog],
                           task_logs: Option[List[WesLog]],
                           outputs: Option[JsObject]
                          )

final case class CromwellCallsMetadata(shardIndex: Option[Int],
                                       commandLine: Option[String],
                                       returnCode: Option[Int],
                                       start: Option[String],
                                       end: Option[String],
                                       stdout: Option[String],
                                       stderr: Option[String]
                                      )

object CromwellMetadata {

  def cromwellCallsMetadataEntryToLogEntry(taskName: String, callsMetadata: CromwellCallsMetadata): WesLog = {
    val newTaskName = callsMetadata.shardIndex map {
      case -1 => taskName
      case notMinusOne => s"$taskName.$notMinusOne"
    } getOrElse taskName

    WesLog(
      name = Option(newTaskName),
      cmd = callsMetadata.commandLine.map(c => List(c)),
      start_time = callsMetadata.start,
      end_time = callsMetadata.end,
      stdout = callsMetadata.stdout,
      stderr = callsMetadata.stderr,
      exit_code = callsMetadata.returnCode
    )
  }
}

object WesRunLog {
//  def fromCromwellMetadata(response: Future[MetadataJsonResponse]): WesResponse = {
//
//    response match {
//      case w: SuccessfulMetadataJsonResponse =>
//        val runs = w.responseJson.results.toList.map(x => WesRunLog(x.run_id)
//
//      case e: FailedMetadataJsonResponse =>
//    }
//  }
}