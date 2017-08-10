package centaur.test.submit

import cats.data.Validated._
import cats.implicits._
import centaur.test.ErrorOr
import com.typesafe.config.Config
import configs.Result
import configs.syntax._

case class SubmitResponse(statusCode: Int, message: String)

object SubmitResponse {

  def fromConfig(conf: Config): ErrorOr[Option[SubmitResponse]] = {
    conf.get[Config]("submit") match {
      case Result.Failure(_) => Valid(None)
      case Result.Success(submitConf) =>
        val errorOrStatusCode: ErrorOr[Int] = toErrorOr(submitConf.get[Int]("statusCode"))
        val errorOrMessage: ErrorOr[String] = toErrorOr(submitConf.get[String]("message"))
        errorOrStatusCode |@| errorOrMessage map {
          SubmitResponse(_, _)
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
