package centaur.test.standard

import cats.data.Validated.{Valid, _}
import centaur.test.standard.StandardTestFormat.{WorkflowFailureTest, WorkflowSuccessTest}
import centaur.test.ErrorOr
import com.typesafe.config.Config
import configs.Result.{Failure, Success}
import configs.syntax._

sealed abstract class StandardTestFormat(val name: String) {
  def testSpecString = this match {
    case WorkflowSuccessTest => "successfully run"
    case WorkflowFailureTest => "fail during execution"
  }
}

object StandardTestFormat {

  case object WorkflowSuccessTest extends StandardTestFormat("WorkflowSuccess")
  case object WorkflowFailureTest extends StandardTestFormat("WorkflowFailure")

  def fromConfig(conf: Config): ErrorOr[StandardTestFormat] = {
    conf.get[String]("testFormat") match {
      case Success(f) => StandardTestFormat.fromString(f)
      case Failure(_) => invalidNel("No testFormat string provided")
    }
  }

  def fromString(testFormat: String): ErrorOr[StandardTestFormat] = {
    val formats = List(WorkflowSuccessTest, WorkflowFailureTest)
    formats collectFirst {
      case format if format.name.equalsIgnoreCase(testFormat) => Valid(format)
    } getOrElse invalidNel(s"No such test format: $testFormat")
  }
}
