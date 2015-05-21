package cromwell

import akka.actor.ActorSystem
import akka.testkit.{EventFilter, DefaultTimeout, ImplicitSender, TestKit}
import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}

abstract class CromwellSpec(actorSystem: ActorSystem) extends TestKit(actorSystem) with DefaultTimeout
with ImplicitSender with WordSpecLike with Matchers with BeforeAndAfterAll {
  
  def startingCallsFilter(callNames: String*): EventFilter =
    EventFilter.info(message = "Starting calls: " + callNames.mkString(", "), occurrences = 1)
  
}
