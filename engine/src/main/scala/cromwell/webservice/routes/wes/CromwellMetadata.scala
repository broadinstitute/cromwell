package cromwell.webservice.routes.wes

import spray.json.{JsObject, JsonFormat, JsonParser}

final case class CromwellSubmittedFiles(workflow: Option[String],
                                        workflowType: Option[String],
                                        workflowTypeVersion: Option[String],
                                        options: Option[String],
                                        inputs: Option[String],
                                        labels: Option[String]
)

final case class CromwellCallsMetadata(shardIndex: Option[Int],
                                       commandLine: Option[String],
                                       returnCode: Option[Int],
                                       start: Option[String],
                                       end: Option[String],
                                       stdout: Option[String],
                                       stderr: Option[String]
)

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

    val workflowRequest = WesRunRequest(
      workflow_params = workflowParams,
      workflow_type = submittedFiles.workflowType.getOrElse("None supplied"),
      workflow_type_version = submittedFiles.workflowTypeVersion.getOrElse("None supplied"),
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
      state = WesState.fromStatusString(Option(status)),
      run_log = Option(workflowLogData),
      task_logs = Option(taskLogs),
      outputs = outputs
    )
  }
}

object CromwellMetadata {
  import spray.json.DefaultJsonProtocol._

  implicit val cromwellCallsMetadataFormat: JsonFormat[CromwellCallsMetadata] = jsonFormat7(CromwellCallsMetadata.apply)
  implicit val cromwellSubmittedFilesFormat: JsonFormat[CromwellSubmittedFiles] = jsonFormat6(
    CromwellSubmittedFiles.apply
  )
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
