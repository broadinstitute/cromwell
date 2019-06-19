package centaur.test.standard

import centaur.test.markers.CallMarker
import centaur.test.standard.CentaurTestFormat._
import com.typesafe.config.Config
import common.Checked
import configs.Result.{Failure, Success}
import configs.syntax._

sealed abstract class CentaurTestFormat(val name: String) {
  val lowerCaseName = name.toLowerCase
  
  def testSpecString: String = this match {
    case WorkflowSuccessTest => "successfully run"
    case WorkflowFailureTest => "fail during execution"
    case RunTwiceExpectingCallCachingTest => "call cache the second run of"
    case RunThriceExpectingCallCachingTest => "call cache the third run of"
    case RunTwiceExpectingNoCallCachingTest => "NOT call cache the second run of"
    case RunFailingTwiceExpectingNoCallCachingTest => "Fail the first run and NOT call cache the second run of"
    case SubmitFailureTest => "fail to submit"
    case InstantAbort => "abort a workflow immediately after submission"
    case _: PapiUpgradeTest => "make sure a PAPI v1 to v2 upgrade preserves call caching when the `name-for-call-caching-purposes` attribute is used"
    case _: CromwellRestartWithRecover => "survive a Cromwell restart and recover jobs"
    case _: CromwellRestartWithoutRecover => "survive a Cromwell restart"
    case _: ScheduledAbort => "abort a workflow mid run"
    case _: ScheduledAbortWithRestart => "abort a workflow mid run and restart immediately"
    case _: WorkflowFailureRestartWithRecover => "survive a Cromwell restart when a workflow was failing and recover jobs"
    case _: WorkflowFailureRestartWithoutRecover => "survive a Cromwell restart when a workflow was failing"
    case other => s"unrecognized format $other"
  }

  def isParallel: Boolean = true
}

object CentaurTestFormat {
  import common.validation.Checked._

  sealed trait SequentialTestFormat extends CentaurTestFormat {
    override def isParallel: Boolean = false
  }
 
  sealed trait RestartFormat extends SequentialTestFormat
  sealed trait WithCallMarker { this: CentaurTestFormat => val build: CallMarker => CentaurTestFormat }
  
  case object WorkflowSuccessTest extends CentaurTestFormat("WorkflowSuccess")
  case object WorkflowFailureTest extends CentaurTestFormat("WorkflowFailure")
  case object RunTwiceExpectingCallCachingTest extends CentaurTestFormat("RunTwiceExpectingCallCaching")
  case object RunThriceExpectingCallCachingTest extends CentaurTestFormat(name = "RunThriceExpectingCallCaching")
  case object RunTwiceExpectingNoCallCachingTest extends CentaurTestFormat("RunTwiceExpectingNoCallCaching")
  case object RunFailingTwiceExpectingNoCallCachingTest extends CentaurTestFormat("RunFailingTwiceExpectingNoCallCaching")
  case object SubmitFailureTest extends CentaurTestFormat("SubmitFailure")
  case object InstantAbort extends CentaurTestFormat("InstantAbort") with SequentialTestFormat

  object CromwellRestartWithRecover extends CentaurTestFormat("CromwellRestartWithRecover") with WithCallMarker {
    val build = CromwellRestartWithRecover.apply _
  }
  case class CromwellRestartWithRecover(callMarker: CallMarker) extends CentaurTestFormat(CromwellRestartWithRecover.name) with RestartFormat
  
  object CromwellRestartWithoutRecover extends CentaurTestFormat("CromwellRestartWithoutRecover") with WithCallMarker {
    val build = CromwellRestartWithoutRecover.apply _
  }
  case class CromwellRestartWithoutRecover(callMarker: CallMarker) extends CentaurTestFormat(CromwellRestartWithoutRecover.name) with RestartFormat

  object ScheduledAbort extends CentaurTestFormat("ScheduledAbort") with WithCallMarker {
    val build = ScheduledAbort.apply _
  }
  case class ScheduledAbort(callMarker: CallMarker) extends CentaurTestFormat(ScheduledAbort.name) with SequentialTestFormat

  object ScheduledAbortWithRestart extends CentaurTestFormat("ScheduledAbortWithRestart") with WithCallMarker {
    val build = ScheduledAbortWithRestart.apply _
  }
  case class ScheduledAbortWithRestart(callMarker: CallMarker) extends CentaurTestFormat(ScheduledAbortWithRestart.name) with RestartFormat

  object WorkflowFailureRestartWithRecover extends CentaurTestFormat("WorkflowFailureRestartWithRecover") with WithCallMarker {
    val build = WorkflowFailureRestartWithRecover.apply _
  }
  case class WorkflowFailureRestartWithRecover(callMarker: CallMarker) extends CentaurTestFormat(WorkflowFailureRestartWithRecover.name) with RestartFormat

  object WorkflowFailureRestartWithoutRecover extends CentaurTestFormat("WorkflowFailureRestartWithoutRecover") with WithCallMarker {
    val build = WorkflowFailureRestartWithoutRecover.apply _
  }
  case class WorkflowFailureRestartWithoutRecover(callMarker: CallMarker) extends CentaurTestFormat(WorkflowFailureRestartWithoutRecover.name) with RestartFormat

  object PapiUpgradeTest extends CentaurTestFormat("PapiUpgrade") with WithCallMarker {
    val build = PapiUpgradeTest.apply _
  }
  case class PapiUpgradeTest(callMarker: CallMarker) extends CentaurTestFormat(PapiUpgradeTest.name) with RestartFormat

  def fromConfig(conf: Config): Checked[CentaurTestFormat] = {
    
    CallMarker.fromConfig(conf).toEither flatMap { callMarker =>
      conf.get[String]("testFormat") match {
        case Success(f) => CentaurTestFormat.fromString(f, callMarker)
        case Failure(_) => "No testFormat string provided".invalidNelCheck[CentaurTestFormat]
      }
    }
  }

  private def fromString(testFormat: String, callMarker: Option[CallMarker]): Checked[CentaurTestFormat] = {
    def withCallMarker(name: String, constructor: CallMarker => CentaurTestFormat) = callMarker match {
      case Some(marker) => constructor(marker).validNelCheck
      case None => s"$name needs a callMarker to know on which call to trigger the restart".invalidNelCheck[CentaurTestFormat]
    }
    
    List(
      WorkflowSuccessTest,
      WorkflowFailureTest,
      RunTwiceExpectingCallCachingTest,
      RunThriceExpectingCallCachingTest,
      RunTwiceExpectingNoCallCachingTest,
      RunFailingTwiceExpectingNoCallCachingTest,
      CromwellRestartWithRecover,
      CromwellRestartWithoutRecover,
      SubmitFailureTest,
      ScheduledAbortWithRestart,
      InstantAbort,
      ScheduledAbort,
      WorkflowFailureRestartWithRecover,
      WorkflowFailureRestartWithoutRecover,
      PapiUpgradeTest
    ).collectFirst({
      case format: WithCallMarker if format.name.equalsIgnoreCase(testFormat) => withCallMarker(format.name, format.build)
      case format if format.name.equalsIgnoreCase(testFormat) => format.validNelCheck
    }).getOrElse(s"No such test format: $testFormat".invalidNelCheck[CentaurTestFormat])
  }
}
