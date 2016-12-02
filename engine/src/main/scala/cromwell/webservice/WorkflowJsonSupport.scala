package cromwell.webservice

import java.nio.file.Paths
import java.time.OffsetDateTime

import cromwell.core._
import cromwell.engine._
import cromwell.services.metadata.MetadataService
import MetadataService.{WorkflowQueryResponse, WorkflowQueryResult}
import cromwell.util.JsonFormatting.WdlValueJsonFormatter
import WdlValueJsonFormatter._
import better.files.File
import spray.json.{DefaultJsonProtocol, JsString, JsValue, RootJsonFormat}

object WorkflowJsonSupport extends DefaultJsonProtocol {
  implicit val workflowStatusResponseProtocol = jsonFormat2(WorkflowStatusResponse)
  implicit val workflowAbortResponseProtocol = jsonFormat2(WorkflowAbortResponse)
  implicit val workflowSubmitResponseProtocol = jsonFormat2(WorkflowSubmitResponse)
  implicit val workflowOutputResponseProtocol = jsonFormat2(WorkflowOutputResponse)
  implicit val callOutputResponseProtocol = jsonFormat3(CallOutputResponse)
  implicit val engineStatsProtocol = jsonFormat2(EngineStatsActor.EngineStats)
  implicit val callAttempt = jsonFormat2(CallAttempt)
  implicit val workflowSourceData = jsonFormat3(WorkflowSourceFilesWithoutImports)

  implicit object fileJsonFormat extends RootJsonFormat[File] {
    override def write(obj: File) = JsString(obj.path.toAbsolutePath.toString)
    override def read(json: JsValue): File = json match {
      case JsString(str) => Paths.get(str)
      case unknown => throw new NotImplementedError(s"Cannot parse $unknown to a File")
    }
  }

  implicit val workflowSourceDataWithImports = jsonFormat4(WorkflowSourceFilesWithDependenciesZip)
  implicit val errorResponse = jsonFormat3(FailureResponse)
  implicit val successResponse = jsonFormat3(SuccessResponse)

  implicit object DateJsonFormat extends RootJsonFormat[OffsetDateTime] {
    override def write(obj: OffsetDateTime) = JsString(obj.toString)

    override def read(json: JsValue): OffsetDateTime = json match {
      case JsString(str) => OffsetDateTime.parse(str)
      case unknown => throw new NotImplementedError(s"Cannot parse $unknown to a DateTime")
    }
  }

  implicit val unqualifiedFailureEventEntry = jsonFormat2(FailureEventEntry)
  implicit val workflowQueryResult = jsonFormat5(WorkflowQueryResult)
  implicit val workflowQueryResponse = jsonFormat1(WorkflowQueryResponse)
}

