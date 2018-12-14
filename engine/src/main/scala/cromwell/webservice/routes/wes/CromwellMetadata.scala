package cromwell.webservice.routes.wes

import cromwell.core.WorkflowState
import spray.json.{DefaultJsonProtocol, JsObject, JsonFormat, JsonParser}

/**
  * As described in WesRouteSupport the least bad solution to reshaping the outgoing metadata JSON at the moment
  * is to request Cromwell's metadata JSON, parsing it back into Scala objects, and then outputting back into
  * the correct WES JSON.
  *
  * This class is regretably required at the moment as we don't currently have a Scala based object modeling of Cromwell's
  * metadata response anywhere in the code base. This is not a full representation, just the elements actually
  * needed for WES. Perhaps some day in the future we'll have a proper client library, but hopefully by the time
  * that day arrives this whole json->scala->json scheme won't be necessary anyways
  */
final case class CromwellMetadata(workflowName: Option[String],
                                  id: String,
                                  status: String,
                                  start: Option[String],
                                  end: Option[String],
                                  submittedFiles: CromwellSubmittedFiles,
                                  outputs: Option[JsObject],
                                  actualWorkflowLanguage: Option[String],
                                  actualWorkflowLanguageVersion: Option[String],
                                  calls: Option[Map[String, Seq[CromwellCallsMetadata]]]
                                 ) {
  import CromwellMetadata._

  def wesResponseRunLog: WesResponseRunLog = {
    val workflowParams = submittedFiles.inputs.map(JsonParser(_).asJsObject)
    val workflowTags = submittedFiles.labels.map(JsonParser(_).asJsObject)
    val workflowEngineParams = submittedFiles.options.map(JsonParser(_).asJsObject)

    val workflowRequest = WesRunRequest(workflow_params = workflowParams,
      workflow_type = actualWorkflowLanguage,
      workflow_type_version = actualWorkflowLanguageVersion,
      tags = workflowTags,
      workflow_engine_parameters = workflowEngineParams,
      workflow_url = submittedFiles.workflowUrl
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

    WesResponseRunLog(
      run_id = id,
      request = workflowRequest,
      state = WesState.fromCromwellStatus(WorkflowState.withName(status)),
      run_log = Option(workflowLogData),
      task_logs = Option(taskLogs),
      outputs = outputs
    )
  }
}

final case class CromwellCallsMetadata(shardIndex: Option[Int],
                                       commandLine: Option[String],
                                       returnCode: Option[Int],
                                       start: Option[String],
                                       end: Option[String],
                                       stdout: Option[String],
                                       stderr: Option[String]
                                      )

final case class CromwellSubmittedFiles(workflow: Option[String],
                                        options: Option[String],
                                        inputs: Option[String],
                                        labels: Option[String],
                                        workflowUrl: Option[String]
                                       )

object CromwellMetadata {
  import DefaultJsonProtocol._

  implicit val cromwellCallsMetadataFormat: JsonFormat[CromwellCallsMetadata] = jsonFormat7(CromwellCallsMetadata.apply)
  implicit val cromwellSubmittedFilesFormat: JsonFormat[CromwellSubmittedFiles] = jsonFormat5(CromwellSubmittedFiles.apply)
  implicit val cromwellMetadataFormat: JsonFormat[CromwellMetadata] = jsonFormat10(CromwellMetadata.apply)

  def fromJson(json: JsObject): CromwellMetadata = json.convertTo[CromwellMetadata]

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
