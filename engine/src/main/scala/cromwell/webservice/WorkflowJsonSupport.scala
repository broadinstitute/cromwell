package cromwell.webservice

import java.nio.file.Paths
import java.time.OffsetDateTime

import cromwell.core._
import cromwell.engine._
import cromwell.services.metadata.MetadataService
import MetadataService._
import cromwell.util.JsonFormatting.WomValueJsonFormatter
import WomValueJsonFormatter._
import better.files.File
import cromwell.services.healthmonitor.HealthMonitorServiceActor.{StatusCheckResponse, SubsystemStatus}
import cromwell.webservice.CromwellApiService.BackendResponse
import cromwell.webservice.metadata.MetadataBuilderActor.BuiltMetadataResponse
import spray.json.{DefaultJsonProtocol, JsString, JsValue, RootJsonFormat}

object WorkflowJsonSupport extends DefaultJsonProtocol {
  implicit val workflowStatusResponseProtocol = jsonFormat2(WorkflowStatusResponse)
  implicit val workflowAbortResponseProtocol = jsonFormat2(WorkflowAbortResponse)
  implicit val workflowSubmitResponseProtocol = jsonFormat2(WorkflowSubmitResponse)
  implicit val workflowOutputResponseProtocol = jsonFormat2(WorkflowOutputResponse)
  implicit val callOutputResponseProtocol = jsonFormat3(CallOutputResponse)
  implicit val engineStatsProtocol = jsonFormat2(EngineStatsActor.EngineStats)
  implicit val BackendResponseFormat = jsonFormat2(BackendResponse)
  implicit val BuiltStatusResponseFormat = jsonFormat1(BuiltMetadataResponse)
  implicit val callAttempt = jsonFormat2(CallAttempt)
  implicit val workflowSourceData = jsonFormat9(WorkflowSourceFilesWithoutImports)
  implicit val subsystemStatusFormat = jsonFormat2(SubsystemStatus)
  implicit val statusCheckResponseFormat = jsonFormat2(StatusCheckResponse)

  implicit object fileJsonFormat extends RootJsonFormat[File] {
    override def write(obj: File) = JsString(obj.path.toAbsolutePath.toString)
    override def read(json: JsValue): File = json match {
      case JsString(str) => Paths.get(str)
      case unknown => throw new NotImplementedError(s"Cannot parse $unknown to a File")
    }
  }

  implicit val workflowSourceDataWithImports = jsonFormat10(WorkflowSourceFilesWithDependenciesZip)
  implicit val errorResponse = jsonFormat3(FailureResponse)
  implicit val successResponse = jsonFormat3(SuccessResponse)

  implicit object DateJsonFormat extends RootJsonFormat[OffsetDateTime] {
    override def write(obj: OffsetDateTime) = JsString(obj.toString)

    override def read(json: JsValue): OffsetDateTime = json match {
      case JsString(str) => OffsetDateTime.parse(str)
      case unknown => throw new NotImplementedError(s"Cannot parse $unknown to a DateTime")
    }
  }

  implicit val workflowQueryResult = jsonFormat8(WorkflowQueryResult)
  implicit val workflowQueryResponse = jsonFormat2(WorkflowQueryResponse)
}
