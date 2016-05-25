package centaur.test.standard

import cats.data.Validated.{Valid, _}
import centaur.test.standard.StandardTestFormat.{SubmissionFailureTest, WorkflowFailureTest, WorkflowSuccessTest}
import centaur.test.ErrorOr
import com.typesafe.config.Config
import configs.Result.{Failure, Success}
import configs.syntax._

sealed abstract class StandardTestFormat {
  def testSpecString = this match {
    case WorkflowSuccessTest => "successfully run"
    case WorkflowFailureTest => "fail during execution"
    case SubmissionFailureTest => "fail to submit"
  }
}

object StandardTestFormat {
  case object WorkflowSuccessTest extends StandardTestFormat
  case object WorkflowFailureTest extends StandardTestFormat
  case object SubmissionFailureTest extends StandardTestFormat

  def fromConfig(conf: Config): ErrorOr[StandardTestFormat] = {
    conf.get[String]("testFormat") match {
      case Success(f) => StandardTestFormat.fromString(f)
      case Failure(_) => invalidNel("No testFormat string provided")
    }
  }

  def fromString(testFormat: String): ErrorOr[StandardTestFormat] = {
    testFormat.toLowerCase match {
      case "workflowsuccess" => Valid(WorkflowSuccessTest)
      case "workflowfailure" => Valid(WorkflowFailureTest)
      case "submissionfailure" => Valid(SubmissionFailureTest)
      case bad => invalidNel(s"No such test format: $bad")
    }
  }
}
