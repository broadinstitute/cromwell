package centaur.test.standard

import cats.data.Validated.{Valid, _}
import centaur.test.{ErrorOr, TestFormulas}
import com.typesafe.config.Config
import lenthall.config.ScalaConfig._

sealed abstract class StandardTestFormat {
  def testSpecString = this match {
    case WorkflowSuccessTest => "successfully run"
    case WorkflowFailureTest => "fail during execution"
    case SubmissionFailureTest => "fail to submit"
  }

  def testFunction = this match {
    case WorkflowSuccessTest => TestFormulas.runSuccessfulWorkflow _
    case WorkflowFailureTest => TestFormulas.runFailingWorkflow _
    case SubmissionFailureTest => TestFormulas.runSubmissionFailureWorkflow _
  }
}

case object WorkflowSuccessTest extends StandardTestFormat
case object WorkflowFailureTest extends StandardTestFormat
case object SubmissionFailureTest extends StandardTestFormat

object StandardTestFormat {
  def fromConfig(conf: Config): ErrorOr[StandardTestFormat] = {
    conf.getStringOption("testFormat") match {
      case Some(f) => StandardTestFormat.fromString(f)
      case None => invalidNel("No testFormat string provided")
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
