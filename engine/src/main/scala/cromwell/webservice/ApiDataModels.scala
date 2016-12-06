package cromwell.webservice

import lenthall.exception.MessageAggregation
import spray.json._
import wdl4s.FullyQualifiedName
import wdl4s.values.WdlValue


case class WorkflowStatusResponse(id: String, status: String)

case class WorkflowSubmitResponse(id: String, status: String)

case class WorkflowOutputResponse(id: String, outputs: Map[FullyQualifiedName, WdlValue])

case class WorkflowAbortResponse(id: String, status: String)

case class CallOutputResponse(id: String, callFqn: String, outputs: Map[FullyQualifiedName, WdlValue])

case class WorkflowMetadataQueryParameters(outputs: Boolean = true, timings: Boolean = true)

object APIResponse {

  private def constructFailureResponse(status: String, ex: Throwable) ={
    ex match {
      case exceptionWithErrors: MessageAggregation =>
        FailureResponse(status, exceptionWithErrors.getMessage,
          Option(JsArray(exceptionWithErrors.errorMessages.toList.map(JsString(_)).toVector)))
      case e: Throwable => FailureResponse(status, e.getMessage, None)
    }
  }

  /** When the data submitted in the request is incorrect. */
  def fail(ex: Throwable) = constructFailureResponse("fail", ex)

  /** When an exception is thrown while processing the request. */
  def error(ex: Throwable) = constructFailureResponse("error", ex)

  /** When a request completes successfully. */
  def success(message: String, data: Option[JsValue] = None) = SuccessResponse("success", message, data)
}

case class SuccessResponse(status: String, message: String, data: Option[JsValue])
case class FailureResponse(status: String, message: String, errors: Option[JsValue])
