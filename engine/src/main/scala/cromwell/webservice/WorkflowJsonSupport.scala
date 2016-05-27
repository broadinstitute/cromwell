package cromwell.webservice

import java.time.OffsetDateTime
import java.time.format.DateTimeFormatter

import cromwell.backend.ExecutionEventEntry
import cromwell.engine._
import cromwell.engine.backend.{CallLogs, OldStyleCallMetadata, WorkflowQueryResult}
import cromwell.engine.db.ExecutionDatabaseKey
import cromwell.engine.workflow.CallCacheData
import cromwell.webservice.WdlFileJsonFormatter._
import cromwell.webservice.WdlValueJsonFormatter._
import spray.json.{DefaultJsonProtocol, JsString, JsValue, RootJsonFormat}

object WorkflowJsonSupport extends DefaultJsonProtocol {
  implicit val workflowStatusResponseProtocol = jsonFormat2(WorkflowStatusResponse)
  implicit val workflowAbortResponseProtocol = jsonFormat2(WorkflowAbortResponse)
  implicit val workflowSubmitResponseProtocol = jsonFormat2(WorkflowSubmitResponse)
  implicit val workflowOutputResponseProtocol = jsonFormat2(WorkflowOutputResponse)
  implicit val callOutputResponseProtocol = jsonFormat3(CallOutputResponse)
  implicit val callLogsResponseProtocol = jsonFormat3(CallLogs)
  implicit val callAttempt = jsonFormat2(CallAttempt)
  implicit val callStdoutStderrResponse = jsonFormat2(CallStdoutStderrResponse)
  implicit val workflowSourceData = jsonFormat3(WorkflowSourceFiles)
  implicit val errorResponse = jsonFormat3(FailureResponse)
  implicit val successResponse = jsonFormat3(SuccessResponse)

  implicit object DateJsonFormat extends RootJsonFormat[OffsetDateTime] {
    override def write(obj: OffsetDateTime) = JsString(obj.format(DateTimeFormatter.ISO_OFFSET_DATE_TIME))

    override def read(json: JsValue): OffsetDateTime = json match {
      case JsString(str) => OffsetDateTime.parse(str)
      case unknown => throw new NotImplementedError(s"Cannot parse $unknown to a DateTime")
    }
  }
  implicit val executionDatabaseKeyValue = jsonFormat3(ExecutionDatabaseKey)
  implicit val unqualifiedFailureEventEntry = jsonFormat2(FailureEventEntry)
  implicit val qualifiedFailureEventEntry = jsonFormat4(QualifiedFailureEventEntry)
  implicit val executionEventProtocol = jsonFormat3(ExecutionEventEntry)
  implicit val callCacheHitProtocol = jsonFormat3(CallCacheData)
  implicit val callMetadataProtocol = jsonFormat19(OldStyleCallMetadata)
  implicit val workflowMetadataResponse = jsonFormat10(WorkflowMetadataResponse)
  implicit val workflowFailuresResponse = jsonFormat4(WorkflowFailuresResponse)
  implicit val workflowQueryResult = jsonFormat5(WorkflowQueryResult)
  implicit val workflowQueryResponse = jsonFormat1(WorkflowQueryResponse)
  implicit val callCachingResponse = jsonFormat1(CallCachingResponse)
}

