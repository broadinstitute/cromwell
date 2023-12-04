package cromwell.engine.workflow.lifecycle.execution.callcaching

import cats.data.Validated._
import cats.syntax.apply._
import cats.syntax.validated._
import cromwell.core.WorkflowId
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheDiffQueryParameter.CallCacheDiffQueryCall
import common.validation.ErrorOr.{ErrorOr, ShortCircuitingFlatMap}

import scala.util.{Failure, Success, Try}

object CallCacheDiffQueryParameter {
  case class CallCacheDiffQueryCall(workflowId: WorkflowId, callFqn: String, jobIndex: Option[Int])

  private def missingWorkflowError(attribute: String) = s"missing $attribute query parameter".invalidNel

  def fromParameters(parameters: Seq[(String, String)]): ErrorOr[CallCacheDiffQueryParameter] = {
    def extractIndex(parameter: String): ErrorOr[Option[Int]] =
      parameters.find(_._1 == parameter) match {
        case Some((_, value)) =>
          Try(value.trim.toInt) match {
            case Success(index) => Option(index).validNel
            case Failure(f) => f.getMessage.invalidNel
          }
        case None => None.validNel
      }

    def extractAttribute(parameter: String): ErrorOr[String] =
      parameters.find(_._1 == parameter) match {
        case Some((_, value)) => value.validNel
        case None => missingWorkflowError(parameter)
      }

    def validateWorkflowId(parameter: String): ErrorOr[WorkflowId] = for {
      workflowIdString <- extractAttribute(parameter)
      workflowId <- fromTry(Try(WorkflowId.fromString(workflowIdString.trim)))
        .leftMap(_.getMessage)
        .toValidatedNel[String, WorkflowId]
    } yield workflowId

    val workflowAValidation = validateWorkflowId("workflowA")
    val workflowBValidation = validateWorkflowId("workflowB")

    val callAValidation: ErrorOr[String] = extractAttribute("callA").map(_.trim)
    val callBValidation: ErrorOr[String] = extractAttribute("callB").map(_.trim)

    val indexAValidation: ErrorOr[Option[Int]] = extractIndex("indexA")
    val indexBValidation: ErrorOr[Option[Int]] = extractIndex("indexB")

    (workflowAValidation,
     callAValidation,
     indexAValidation,
     workflowBValidation,
     callBValidation,
     indexBValidation
    ) mapN { (workflowA, callA, indexA, workflowB, callB, indexB) =>
      CallCacheDiffQueryParameter(
        CallCacheDiffQueryCall(workflowA, callA, indexA),
        CallCacheDiffQueryCall(workflowB, callB, indexB)
      )
    }
  }
}

case class CallCacheDiffQueryParameter(callA: CallCacheDiffQueryCall, callB: CallCacheDiffQueryCall)
