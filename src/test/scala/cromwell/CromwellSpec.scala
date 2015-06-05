package cromwell

import akka.actor.ActorSystem
import akka.testkit.{DefaultTimeout, EventFilter, ImplicitSender, TestKit}
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}

abstract class CromwellSpec(actorSystem: ActorSystem) extends TestKit(actorSystem) with DefaultTimeout
with ImplicitSender with WordSpecLike with Matchers with BeforeAndAfterAll with ScalaFutures {

  def startingCallsFilter(callNames: String*): EventFilter =
    EventFilter.info(message = s"Starting calls: ${callNames.mkString(", ")}", occurrences = 1)

  def waitForHandledMessage[T](named: String)(block: => T): T =
    EventFilter.debug(start = s"received handled message $named", occurrences = 1).intercept {
      block
    }
}
