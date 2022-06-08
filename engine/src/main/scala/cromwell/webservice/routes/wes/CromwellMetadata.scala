package cromwell.webservice.routes.wes

import WesState._
import spray.json.{JsObject, JsonFormat, JsonParser}

object CromwellMetadata {
  import spray.json.DefaultJsonProtocol._

  implicit val cromwellMetadataFormat: JsonFormat[CromwellMetadata] = jsonFormat8(CromwellMetadata.apply)

  def fromJson(json: String): CromwellMetadata = {
    val jsonAst = JsonParser(json)
    jsonAst.convertTo[CromwellMetadata]
  }

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

final case class CromwellMetadata(workflowName: Option[String],
                                  id: String,
                                  status: String,
                                  start: Option[String],
                                  end: Option[String],
                                  submittedFiles: CromwellSubmittedFiles,
                                  outputs: Option[JsObject],
                                  calls: Option[Map[String, Seq[CromwellCallsMetadata]]]
                                 ) {
  import CromwellMetadata._

  def wesRunLog: WesRunLog = {
    val workflowParams = submittedFiles.inputs.map(JsonParser(_).asJsObject)
    val workflowTags = submittedFiles.labels.map(JsonParser(_).asJsObject)
    val workflowEngineParams = submittedFiles.options.map(JsonParser(_).asJsObject)

    val workflowRequest = WesRunRequest(workflow_params = workflowParams,
      workflow_type = submittedFiles.workflowType.getOrElse("Unable to find workflow type"),
      workflow_type_version = submittedFiles.workflowTypeVersion.getOrElse("Unable to find workflow version"),
      tags = workflowTags,
      workflow_engine_parameters = workflowEngineParams,
      workflow_url = None
    )

    val workflowLogData = WesLog(name = workflowName,
      cmd = None,
      start_time = start,
      end_time = end,
      stdout = None,
      stderr = None,
      exit_code = None
    )

    val taskLogs = for {
      callsArray <- calls.toList
      (taskName, metadataEntries) <- callsArray
      metadataEntry <- metadataEntries
      logEntry = cromwellCallsMetadataEntryToLogEntry(taskName, metadataEntry)
    } yield logEntry

    WesRunLog(
      run_id = id,
      request = workflowRequest,
      state = WesState.fromCromwellStatus(status),
      run_log = Option(workflowLogData),
      task_logs = Option(taskLogs),
      outputs = outputs
    )
  }
}
