package cromwell.engine.workflow.lifecycle.execution.callcaching

import cats.data.{NonEmptyList, Validated}
import cats.implicits._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheDiffQueryParameter.CallCacheDiffQueryCall

import scala.util.{Failure, Success, Try}

object CallCacheDiffQueryParameter {
  case class CallCacheDiffQueryCall(workflowId: String, callFqn: String, jobIndex: Option[Int])

  private def missingWorkflowError(attribute: String) = s"missing $attribute query parameter".invalidNel

  private def extractAttribute(parameters: Seq[(String, String)], parameter: String): Validated[NonEmptyList[String], String] = {
    parameters.find(_._1 == parameter) match {
      case Some((_, value)) => value.validNel
      case None => missingWorkflowError(parameter)
    }
  }

  private def extractIndex(parameters: Seq[(String, String)], parameter: String): Validated[NonEmptyList[String], Option[Int]] = {
    parameters.find(_._1 == parameter) match {
      case Some((_, value)) => Try(value.toInt) match {
        case Success(index) => Option(index).validNel
        case Failure(f) => f.getMessage.invalidNel
      }
      case None => None.validNel
    }
  }

  def fromParameters(parameters: Seq[(String, String)]): Validated[NonEmptyList[String], CallCacheDiffQueryParameter] = {
    val workflowAValidation = extractAttribute(parameters, "workflowA")
    val workflowBValidation = extractAttribute(parameters, "workflowB")

    val callAValidation = extractAttribute(parameters, "callA")
    val callBValidation = extractAttribute(parameters, "callB")

    val indexAValidation = extractIndex(parameters, "indexA")
    val indexBValidation = extractIndex(parameters, "indexB")

    (workflowAValidation |@| callAValidation |@| indexAValidation |@|
      workflowBValidation |@| callBValidation |@| indexBValidation) map {
      case ((workflowA, callA, indexA, workflowB, callB, indexB)) =>
        CallCacheDiffQueryParameter(
          CallCacheDiffQueryCall(workflowA, callA, indexA),
          CallCacheDiffQueryCall(workflowB, callB, indexB)
        )
    }
  }
}

case class CallCacheDiffQueryParameter(callA: CallCacheDiffQueryCall, callB: CallCacheDiffQueryCall)
