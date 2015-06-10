package cromwell

import akka.actor.ActorSystem
import akka.testkit.{DefaultTimeout, EventFilter, ImplicitSender, TestKit}
import com.typesafe.config.ConfigFactory
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}

object CromwellTestkitSpec {
  val akkaConfigString =
    """
      |akka {
      |  loggers = ["akka.event.slf4j.Slf4jLogger", "akka.testkit.TestEventListener"]
      |  loglevel = "DEBUG"
      |  actor {
      |    debug {
      |       receive = on
      |    }
      |  }
      |}
    """.stripMargin
}

abstract class CromwellTestkitSpec(name: String)
  extends TestKit(ActorSystem(name, ConfigFactory.parseString(CromwellTestkitSpec.akkaConfigString)))
  with DefaultTimeout with ImplicitSender with WordSpecLike with Matchers with BeforeAndAfterAll with ScalaFutures {

  def startingCallsFilter(callNames: String*): EventFilter =
    EventFilter.info(pattern = s"starting ${callNames.mkString(", ")}", occurrences = 1)

  def waitForHandledMessage[T](named: String)(block: => T): T = {
    waitForHandledMessagePattern(s"^received handled message $named")(block)
  }

  def waitForHandledMessagePattern[T](pattern: String)(block: => T): T = {
    EventFilter.debug(pattern=pattern, occurrences = 1).intercept {
      block
    }
  }
}
