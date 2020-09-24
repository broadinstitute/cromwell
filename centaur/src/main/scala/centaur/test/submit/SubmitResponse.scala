package centaur.test.submit

import cats.data.Validated._
import cats.syntax.all._
import com.typesafe.config.Config
import common.validation.ErrorOr.ErrorOr
import configs.Result
import configs.syntax._
import cromwell.api.model.SubmittedWorkflow

/**
  * Wraps a response from the cromwell client, either a submitted workflow or an HTTP status/message.
  */
sealed trait SubmitResponse

object SubmitResponse {
  def apply(submittedWorkflow: SubmittedWorkflow): SubmitResponse = {
    SubmitWorkflowResponse(submittedWorkflow)
  }

  def apply(statusCode: Int, message: String): SubmitResponse = {
    SubmitHttpResponse(statusCode, message)
  }
}

case class SubmitWorkflowResponse(submittedWorkflow: SubmittedWorkflow) extends SubmitResponse

case class SubmitHttpResponse(statusCode: Int, message: String) extends SubmitResponse

object SubmitHttpResponse {

  def fromConfig(conf: Config): ErrorOr[Option[SubmitHttpResponse]] = {
    conf.get[Config]("submit") match {
      case Result.Failure(_) => Valid(None)
      case Result.Success(submitConf) =>
        val errorOrStatusCode: ErrorOr[Int] = toErrorOr(submitConf.get[Int]("statusCode"))
        val errorOrMessage: ErrorOr[String] = toErrorOr(submitConf.get[String]("message"))
        (errorOrStatusCode, errorOrMessage) mapN {
          SubmitHttpResponse(_, _)
        } map {
          Option(_)
        }
    }
  }

  private def toErrorOr[A](result: Result[A]): ErrorOr[A] = {
    result match {
      case Result.Success(value) => Valid(value)
      case Result.Failure(error) =>
        error.messages
          .toList
          .toNel
          .getOrElse(throw new RuntimeException("Paranoia... error.messages is a Nel exposed as a Seq."))
          .invalid
    }
  }
}
