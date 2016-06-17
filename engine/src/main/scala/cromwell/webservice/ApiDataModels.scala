package cromwell.webservice

import cromwell.core.WorkflowId
import cromwell.webservice.PerRequest.RequestComplete
import spray.http.StatusCodes
import spray.json._
import wdl4s.values.WdlValue
import wdl4s.{FullyQualifiedName, ThrowableWithErrors}

import scala.language.postfixOps

case class WorkflowStatusResponse(id: String, status: String)

case class WorkflowSubmitResponse(id: String, status: String)

case class WorkflowOutputResponse(id: String, outputs: Map[FullyQualifiedName, WdlValue])

case class WorkflowAbortResponse(id: String, status: String)

case class CallOutputResponse(id: String, callFqn: String, outputs: Map[FullyQualifiedName, WdlValue])

case class WorkflowMetadataQueryParameters(outputs: Boolean = true, timings: Boolean = true)

object APIResponse {
  import WorkflowJsonSupport._
  import spray.httpx.SprayJsonSupport._

  private def constructFailureResponse(status: String, ex: Throwable) ={
    ex match {
      case cex: ThrowableWithErrors => FailureResponse(status, cex.message, Option(JsArray(cex.errors.list.map(JsString(_)).toVector)))
      case e: Throwable => FailureResponse(status, e.getMessage, None)
    }
  }

  /** When the data submitted in the request is incorrect. */
  def fail(ex: Throwable) = constructFailureResponse("fail", ex)

  /** When an exception is thrown while processing the request. */
  def error(ex: Throwable) = constructFailureResponse("error", ex)

  /** When a request completes successfully. */
  def success(message: String, data: Option[JsValue] = None) = SuccessResponse("success", message, data)

  def workflowNotFound(id: WorkflowId) = RequestComplete(StatusCodes.NotFound, APIResponse.error(new Exception(s"Workflow '$id' not found.")))
}

case class SuccessResponse(status: String, message: String, data: Option[JsValue])
case class FailureResponse(status: String, message: String, errors: Option[JsValue])
