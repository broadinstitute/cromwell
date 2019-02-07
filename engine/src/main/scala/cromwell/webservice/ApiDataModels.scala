package cromwell.webservice

import common.exception.MessageAggregation
import org.apache.commons.lang3.exception.ExceptionUtils
import spray.json._
import wdl.draft2.model.FullyQualifiedName
import wom.values.WomValue


case class WorkflowStatusResponse(id: String, status: String)

case class WorkflowSubmitResponse(id: String, status: String)

case class WorkflowOutputResponse(id: String, outputs: Map[FullyQualifiedName, WomValue])

case class WorkflowAbortResponse(id: String, status: String)

case class CallOutputResponse(id: String, callFqn: String, outputs: Map[FullyQualifiedName, WomValue])

case class WorkflowMetadataQueryParameters(outputs: Boolean = true, timings: Boolean = true)

object APIResponse {

  private def constructFailureResponse(status: String, ex: Throwable) = {
    ex match {
      case exceptionWithErrors: MessageAggregation =>
        FailureResponse(status, exceptionWithErrors.getMessage,
          Option(JsArray(exceptionWithErrors.errorMessages.toList.map(JsString(_)).toVector)))
      case e: Throwable => FailureResponse(status, e.getMessage, Option(e.getCause).map(c => JsArray(JsString(ExceptionUtils.getMessage(c)))))
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
