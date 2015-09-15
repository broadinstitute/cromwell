package cromwell.webservice

import cromwell.binding.values.WdlFileJsonFormatter._
import cromwell.binding.values.WdlValueJsonFormatter._
import cromwell.engine.backend.{CallMetadata, StdoutStderr}
import org.joda.time.format.ISODateTimeFormat
import spray.json.{JsValue, JsString, RootJsonFormat, DefaultJsonProtocol}
import org.joda.time.DateTime


object WorkflowJsonSupport extends DefaultJsonProtocol {
  implicit val workflowValidationResponseProtocol = jsonFormat2(WorkflowValidateResponse)
  implicit val workflowStatusResponseProtocol = jsonFormat2(WorkflowStatusResponse)
  implicit val workflowAbortResponseProtocol = jsonFormat2(WorkflowAbortResponse)
  implicit val workflowSubmitResponseProtocol = jsonFormat2(WorkflowSubmitResponse)
  implicit val workflowOutputResponseProtocol = jsonFormat2(WorkflowOutputResponse)
  implicit val callOutputResponseProtocol = jsonFormat3(CallOutputResponse)
  implicit val callLogsResponseProtocol = jsonFormat2(StdoutStderr)
  implicit val callStdoutStderrResponse = jsonFormat2(CallStdoutStderrResponse)

  implicit object DateJsonFormat extends RootJsonFormat[DateTime] {
    private val parserISO = ISODateTimeFormat.dateTime()
    override def write(obj: DateTime) = JsString(parserISO.print(obj))
    override def read(json: JsValue) : DateTime = ???
  }
  implicit val callMetadataProtocol = jsonFormat7(CallMetadata)
  implicit val callMetadataResponse = jsonFormat8(WorkflowMetadataResponse)
}

