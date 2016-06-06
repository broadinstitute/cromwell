package cromwell.webservice

import java.time.OffsetDateTime

import cromwell.core.WorkflowId
import cromwell.engine.backend.{CallLogs, OldStyleCallMetadata}
import cromwell.engine.{FailureEventEntry, QualifiedFailureEventEntry}
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

case class CallStdoutStderrResponse(id: String, logs: Map[FullyQualifiedName, Seq[CallLogs]])

case class WorkflowMetadataQueryParameters(outputs: Boolean = true, timings: Boolean = true)

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class WorkflowMetadataResponse(id: String,
                                    workflowName: String,
                                    status: String,
                                    submission: OffsetDateTime,
                                    start: Option[OffsetDateTime],
                                    end: Option[OffsetDateTime],
                                    inputs: JsObject,
                                    outputs: Option[Map[String, WdlValue]],
                                    calls: Map[String, Seq[OldStyleCallMetadata]])

final case class CallCachingResponse(updateCount: Int)

case class WorkflowFailuresResponse(id: String,
                                    workflowName: String,
                                    status: String,
                                    failures: Seq[QualifiedFailureEventEntry])

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
