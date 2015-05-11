package cromwell

import akka.actor.ActorSystem
import akka.testkit.{TestKit, ImplicitSender, DefaultTimeout}
import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}

abstract class CromwellSpec(actorSystem: ActorSystem) extends TestKit(actorSystem) with DefaultTimeout
with ImplicitSender with WordSpecLike with Matchers with BeforeAndAfterAll {

}
