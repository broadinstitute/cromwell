package cromwell.webservice

import cromwell.engine.{FailureEventEntry, QualifiedFailureEventEntry}
import cromwell.engine.backend.{CallLogs, CallMetadata, WorkflowQueryResult}
import org.joda.time.DateTime
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

case class WorkflowMetadataResponse(id: String,
                                    workflowName: String,
                                    status: String,
                                    submission: DateTime,
                                    start: Option[DateTime],
                                    end: Option[DateTime],
                                    inputs: JsObject,
                                    outputs: Option[Map[String, WdlValue]],
                                    calls: Map[String, Seq[CallMetadata]],
                                    failures: Option[Seq[FailureEventEntry]])

case class WorkflowQueryResponse(results: Seq[WorkflowQueryResult])

final case class CallCachingResponse(updateCount: Int)

case class WorkflowFailuresResponse(id: String,
                                    workflowName: String,
                                    status: String,
                                    failures: Seq[QualifiedFailureEventEntry])

object APIResponse {
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
}

case class SuccessResponse(status: String, message: String, data: Option[JsValue])
case class FailureResponse(status: String, message: String, errors: Option[JsValue])
