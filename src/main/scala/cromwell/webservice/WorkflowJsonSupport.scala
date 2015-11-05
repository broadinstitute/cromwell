package cromwell.webservice

import cromwell.binding.values.WdlFileJsonFormatter._
import cromwell.binding.values.WdlValueJsonFormatter._
import cromwell.engine.backend.{WorkflowQueryResult, CallMetadata, CallLogs}
import org.joda.time.DateTime
import org.joda.time.format.ISODateTimeFormat
import spray.json.{DefaultJsonProtocol, JsString, JsValue, RootJsonFormat}


object WorkflowJsonSupport extends DefaultJsonProtocol {
  implicit val workflowValidationResponseProtocol = jsonFormat2(WorkflowValidateResponse)
  implicit val workflowStatusResponseProtocol = jsonFormat2(WorkflowStatusResponse)
  implicit val workflowAbortResponseProtocol = jsonFormat2(WorkflowAbortResponse)
  implicit val workflowSubmitResponseProtocol = jsonFormat2(WorkflowSubmitResponse)
  implicit val workflowOutputResponseProtocol = jsonFormat2(WorkflowOutputResponse)
  implicit val callOutputResponseProtocol = jsonFormat3(CallOutputResponse)
  implicit val callLogsResponseProtocol = jsonFormat3(CallLogs)
  implicit val callStdoutStderrResponse = jsonFormat2(CallStdoutStderrResponse)

  implicit object DateJsonFormat extends RootJsonFormat[DateTime] {
    private val parserISO = ISODateTimeFormat.dateTime()
    override def write(obj: DateTime) = JsString(parserISO.print(obj))

    override def read(json: JsValue): DateTime = json match {
      case JsString(str) => parserISO.parseDateTime(str)
      case unknown => throw new NotImplementedError(s"Cannot parse $unknown to a DateTime")
    }
  }
  implicit val callMetadataProtocol = jsonFormat13(CallMetadata)
  implicit val workflowMetadataResponse = jsonFormat8(WorkflowMetadataResponse)
  implicit val workflowQueryResult = jsonFormat5(WorkflowQueryResult)
  implicit val workflowQueryResponse = jsonFormat1(WorkflowQueryResponse)
  implicit val callCachingResponse = jsonFormat1(CallCachingResponse)
}

