package cromwell.webservice

import cromwell.binding.values.WdlValueJsonFormatter._
import spray.http.{HttpEntity, MediaTypes}
import spray.httpx.marshalling.Marshaller
import spray.json.DefaultJsonProtocol
import spray.json._

object WorkflowJsonSupport extends DefaultJsonProtocol {
  implicit val workflowStatusResponseProtocol = jsonFormat2(WorkflowStatusResponse)
  implicit val workflowAbortResponseProtocol = jsonFormat2(WorkflowAbortResponse)
  implicit val workflowSubmitResponseProtocol = jsonFormat2(WorkflowSubmitResponse)
  implicit val workflowOutputResponseProtocol = jsonFormat2(WorkflowOutputResponse)
  implicit val callOutputResponseProtocol = jsonFormat3(CallOutputResponse)
  implicit val callStdoutStderrResponseMarshaller =
    Marshaller.of[CallStdoutStderrResponse](MediaTypes.`application/json`) { (value, contentType, ctx) =>
      val x = value.logs.map {case (fqn, logs) => (fqn, Map("stdout" -> logs.stdout.value, "stderr" -> logs.stderr.value))}
      val map = Map(
        "id" -> value.id,
        "logs" -> x
      )
      ctx.marshalTo(HttpEntity(contentType, map.toJson.prettyPrint))
    }
}
