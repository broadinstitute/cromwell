package cromwell.engine.workflow.lifecycle.execution.callcaching

import cats.data.Validated.{Invalid, Valid}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheDiffQueryParameter.CallCacheDiffQueryCall
import cats.implicits._
import lenthall.exception.AggregatedMessageException

import scala.util.{Failure, Success, Try}

object CallCacheDiffQueryParameter {
  case class CallCacheDiffQueryCall(workflowId: String, callFqn: String, jobIndex: Option[Int])

  private def missingWorkflowError(attribute: String) = s"missing $attribute query parameter".invalidNel

  private def extractAttribute(parameters: Seq[(String, String)], parameter: String) = {
    parameters.find(_._1 == parameter) match {
      case Some((_, value)) => value.validNel
      case None => missingWorkflowError(parameter)
    }
  }

  private def extractIndex(parameters: Seq[(String, String)], parameter: String) = {
    parameters.find(_._1 == parameter) match {
      case Some((_, value)) => Try(value.toInt) match {
        case Success(index) => Option(index).validNel
        case Failure(f) => f.getMessage.invalidNel
      }
      case None => None.validNel
    }
  }

  def fromParameters(parameters: Seq[(String, String)]) = {
    val workflowAValidation = extractAttribute(parameters, "workflowA")
    val workflowBValidation = extractAttribute(parameters, "workflowB")

    val callAValidation = extractAttribute(parameters, "callA")
    val callBValidation = extractAttribute(parameters, "callB")

    val indexAValidation = extractIndex(parameters, "indexA")
    val indexBValidation = extractIndex(parameters, "indexB")

    (workflowAValidation |@| callAValidation |@| indexAValidation |@|
      workflowBValidation |@| callBValidation |@| indexBValidation).tupled match {
      case Valid((workflowA, callA, indexA, workflowB, callB, indexB)) =>
        Success(CallCacheDiffQueryParameter(
          CallCacheDiffQueryCall(workflowA, callA, indexA),
          CallCacheDiffQueryCall(workflowB, callB, indexB)
        ))
      case Invalid(errors) => Failure(AggregatedMessageException("Wrong parameters for call cache diff query", errors.toList))
    }
  }
}

case class CallCacheDiffQueryParameter(callA: CallCacheDiffQueryCall, callB: CallCacheDiffQueryCall)