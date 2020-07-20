package cromwell.webservice

import java.nio.file.Paths
import java.time.OffsetDateTime

import better.files.File
import common.util.TimeUtil._
import common.validation.Validation._
import cromwell.core._
import cromwell.engine._
import cromwell.services.healthmonitor.ProtoHealthMonitorServiceActor.{StatusCheckResponse, SubsystemStatus}
import cromwell.services.metadata.MetadataArchiveStatus
import cromwell.services.metadata.MetadataService._
import cromwell.util.JsonFormatting.WomValueJsonFormatter._
import cromwell.webservice.routes.CromwellApiService.BackendResponse
import spray.json.{DefaultJsonProtocol, JsString, JsValue, JsonFormat, RootJsonFormat}

object WorkflowJsonSupport extends DefaultJsonProtocol {
  implicit val workflowStatusResponseProtocol = jsonFormat2(WorkflowStatusResponse)
  implicit val workflowAbortResponseProtocol = jsonFormat2(WorkflowAbortResponse)
  implicit val workflowSubmitResponseProtocol = jsonFormat2(WorkflowSubmitResponse)
  implicit val workflowOutputResponseProtocol = jsonFormat2(WorkflowOutputResponse)
  implicit val callOutputResponseProtocol = jsonFormat3(CallOutputResponse)
  implicit val engineStatsProtocol = jsonFormat2(EngineStatsActor.EngineStats)
  implicit val BackendResponseFormat = jsonFormat2(BackendResponse)
  implicit val callAttempt = jsonFormat2(CallAttempt)

  implicit val workflowOptionsFormatter: JsonFormat[WorkflowOptions] = new JsonFormat[WorkflowOptions]  {
    override def read(json: JsValue): WorkflowOptions = json match {
      case str: JsString => WorkflowOptions.fromJsonString(str.value).get
      case other => throw new UnsupportedOperationException(s"Cannot use ${other.getClass.getSimpleName} value. Expected a workflow options String")
    }
    override def write(obj: WorkflowOptions): JsValue = JsString(obj.asPrettyJson)
  }

  implicit val workflowSourceData = jsonFormat10(WorkflowSourceFilesWithoutImports)
  implicit val subsystemStatusFormat = jsonFormat2(SubsystemStatus)
  implicit val statusCheckResponseFormat = jsonFormat2(StatusCheckResponse)

  implicit object fileJsonFormat extends RootJsonFormat[File] {
    override def write(obj: File) = JsString(obj.path.toAbsolutePath.toString)
    override def read(json: JsValue): File = json match {
      case JsString(str) => Paths.get(str)
      case unknown => throw new UnsupportedOperationException(s"Cannot parse $unknown to a File")
    }
  }

  implicit val workflowSourceDataWithImports = jsonFormat11(WorkflowSourceFilesWithDependenciesZip)
  implicit val errorResponse = jsonFormat3(FailureResponse)

  // By default the formatter for JsValues prints them out ADT-style.
  // In the case of SuccessResponses, we just want raw JsValues to be included in our output verbatim.
  private implicit val identityJsValueFormatter = new RootJsonFormat[JsValue] {
    override def read(json: JsValue): JsValue = json
    override def write(obj: JsValue): JsValue = obj
  }
  implicit val successResponse = jsonFormat3(SuccessResponse)

  implicit object DateJsonFormat extends RootJsonFormat[OffsetDateTime] {
    override def write(offsetDateTime: OffsetDateTime) = JsString(offsetDateTime.toUtcMilliString)

    override def read(json: JsValue): OffsetDateTime = json match {
      case JsString(str) => OffsetDateTime.parse(str)
      case unknown => throw new UnsupportedOperationException(s"Cannot parse $unknown to a DateTime")
    }
  }

  implicit object MetadataArchiveStatusFormat extends RootJsonFormat[MetadataArchiveStatus] {
    override def write(obj: MetadataArchiveStatus): JsValue = JsString(obj.toString)

    override def read(json: JsValue): MetadataArchiveStatus = json match {
      case JsString(str) => MetadataArchiveStatus.withName(str).toTry.get
      case unknown => throw new UnsupportedOperationException(s"Cannot parse $unknown to a MetadataArchiveStatus")
    }
  }

  implicit val workflowQueryResult = jsonFormat10(WorkflowQueryResult)
  implicit val workflowQueryResponse = jsonFormat2(WorkflowQueryResponse)
}
