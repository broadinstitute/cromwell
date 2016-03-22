package cromwell.core.eventbus

import akka.actor._
import akka.testkit.{DefaultTimeout, ImplicitSender, TestDuration, TestKit}
import cromwell.core.eventbus.TestEntityActor.{Subscribe, Unsubscribe, PublishBye, PublishHello}
import org.scalatest.{BeforeAndAfter, Matchers, WordSpecLike}

import scala.concurrent.duration._
import scala.language.postfixOps

class PubSubMediatorSpec extends TestKit(ActorSystem("PubSubMediatorSpecSystem"))
  with DefaultTimeout with ImplicitSender with WordSpecLike  with BeforeAndAfter with PubSubMediator with Matchers {

  val Timeout = 1.second.dilated
  var testEntityActor: ActorRef = _

  before {
    testEntityActor = system.actorOf(TestEntityActor.props())
  }

  after {
    system.stop(testEntityActor)
  }

  "A PubSubMediator" should {
    "provide the functionality to entities to subscribe to events" in {
      within(Timeout) {
        subscribe[TestTopic](self)
        testEntityActor ! PublishHello
        expectMsg(EventMsg[TestTopic](HelloTopic, "Hello All!"))
      }
    }

    "provide the functionality to entities to unsubscribe from events" in {
      within(Timeout) {
        subscribe[TestTopic](self)
        testEntityActor ! PublishHello
        expectMsg(EventMsg[TestTopic](HelloTopic, "Hello All!"))
        unsubscribe[TestTopic](self)
        testEntityActor ! PublishHello
        expectNoMsg()
      }
    }

    "provide the functionality to entities to publish events" in {
      within(Timeout) {
        testEntityActor ! Subscribe
        subscribe[TestTopic](self)
        publish[TestTopic](EventMsg(AskTopic, "Someone, please respond!"))
        expectMsg(EventMsg(AskTopic, "Someone, please respond!"))
        expectMsg(EventMsg[TestTopic](ReplyTopic, "Yes, sure!"))
      }
    }
  }
}

sealed trait TestTopic
case object HelloTopic extends TestTopic
case object AskTopic extends TestTopic
case object ReplyTopic extends TestTopic

object TestEntityActor {
  sealed trait TestEntityActorMessage
  case object Subscribe extends TestEntityActorMessage
  case object Unsubscribe extends TestEntityActorMessage
  case object PublishHello extends TestEntityActorMessage
  case object PublishBye extends TestEntityActorMessage

  def props(): Props = Props(new TestEntityActor())
}

class TestEntityActor extends Actor with PubSubMediator {
  implicit val actorSystem = context.system

  def receive: Receive = {
    case Subscribe => subscribe[TestTopic](self)
    case Unsubscribe => unsubscribe[TestTopic](self)
    case PublishHello => publish[TestTopic](EventMsg(HelloTopic, "Hello All!"))
    case PublishBye => publish[TestTopic](EventMsg(HelloTopic, "Bye All!"))
    case EventMsg(AskTopic, _) => publish[TestTopic](EventMsg(ReplyTopic, "Yes, sure!"))
  }
}
