package centaur.test.standard

import cats.data.Validated.{Valid, _}
import centaur.test.standard.CentaurTestFormat._
import centaur.test.ErrorOr
import com.typesafe.config.Config
import configs.Result.{Failure, Success}
import configs.syntax._

sealed abstract class CentaurTestFormat(val name: String) {
  def testSpecString: String = this match {
    case WorkflowSuccessTest => "successfully run"
    case WorkflowFailureTest => "fail during execution"
    case RunTwiceExpectingCallCachingTest => "call cache the second run of"
    case RunTwiceExpectingNoCallCachingTest => "NOT call cache the second run of"
    case RunFailingTwiceExpectingNoCallCachingTest => "Fail the first run and NOT call cache the second run of"
    case CromwellRestartWithRecover => "survive a Cromwell restart and recover jobs"
    case CromwellRestartWithoutRecover => "survive a Cromwell restart"
    case SubmitFailureTest => "fail to submit"
  }
}

object CentaurTestFormat {

  case object WorkflowSuccessTest extends CentaurTestFormat("WorkflowSuccess")
  case object WorkflowFailureTest extends CentaurTestFormat("WorkflowFailure")
  case object RunTwiceExpectingCallCachingTest extends CentaurTestFormat("RunTwiceExpectingCallCaching")
  case object RunTwiceExpectingNoCallCachingTest extends CentaurTestFormat("RunTwiceExpectingNoCallCaching")
  case object RunFailingTwiceExpectingNoCallCachingTest extends CentaurTestFormat("RunFailingTwiceExpectingNoCallCaching")
  case object CromwellRestartWithRecover extends CentaurTestFormat("CromwellRestartWithRecover")
  case object CromwellRestartWithoutRecover extends CentaurTestFormat("CromwellRestartWithoutRecover")
  case object SubmitFailureTest extends CentaurTestFormat("SubmitFailure")

  def fromConfig(conf: Config): ErrorOr[CentaurTestFormat] = {
    conf.get[String]("testFormat") match {
      case Success(f) => CentaurTestFormat.fromString(f)
      case Failure(_) => invalidNel("No testFormat string provided")
    }
  }

  def fromString(testFormat: String): ErrorOr[CentaurTestFormat] = {
    val formats = List(
      WorkflowSuccessTest,
      WorkflowFailureTest,
      RunTwiceExpectingCallCachingTest,
      RunTwiceExpectingNoCallCachingTest,
      RunFailingTwiceExpectingNoCallCachingTest,
      CromwellRestartWithRecover,
      CromwellRestartWithoutRecover,
      SubmitFailureTest)
    formats collectFirst {
      case format if format.name.equalsIgnoreCase(testFormat) => Valid(format)
    } getOrElse invalidNel(s"No such test format: $testFormat")
  }
}
