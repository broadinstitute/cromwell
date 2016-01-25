package cromwell.logging

import akka.actor.{Props, ActorRef, ActorSystem}
import akka.testkit.{TestActorRef, TestKit}
import cromwell.pubsub.PubSubMediator
import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}

/**
  * Created by himanshu on 1/21/16.
  */
class BusinessLoggingActorSpec(_system:ActorSystem) extends TestKit(_system) with WordSpecLike with Matchers with BeforeAndAfterAll  {

  var logMessage = "message"
  var subscribed  = false

  trait MockBusinessLogEvent extends BusinessLogEvent {

    override def onBusinessLoggingEvent = (topic: WorkflowEvent, payload: Any) => Some(topic) collect {
      case event@(CallExecutionEvent | WorkflowExecutionEvent) =>
        logMessage = payload.toString
    }
  }

  trait MockPubSubMediator extends PubSubMediator {

    def isSubscribed : Boolean = subscribed

    override def subscribe[T](actorRef: ActorRef)(implicit system:ActorSystem) : Unit= {
      subscribed = true
    }

    override def unsubscribe[T](actorRef: ActorRef,topic:Option[T] = None)(implicit system:ActorSystem) : Unit= {
      subscribed = false
    }
  }


  override def afterAll {
    TestKit.shutdownActorSystem(system)
  }

  def this() = this(ActorSystem("BusinessLoggingActorSpec"))

  "BusinessLogging Actor " must {
    "subscribe a subscriber to " in {
      val actor = TestActorRef(new BusinessLogging() with MockBusinessLogEvent with MockPubSubMediator)
      actor ! SubscribeToLogging
      assert(subscribed == true)
    }
  }

  "BusinessLogging Actor " must {
    "ubsubscribe a subscriber to " in {
      val actor = TestActorRef(new BusinessLogging() with MockBusinessLogEvent with MockPubSubMediator)
      actor ! SubscribeToLogging
      actor ! UnSubscribeToLogging
      assert(subscribed == false)
    }
  }

  "BusinessLogging Actor " must {
    "subscribe a subscriber and pubsubmediator publish message " in {
      val actor = TestActorRef(new BusinessLogging() with MockBusinessLogEvent with PubSubMediator)
      actor ! SubscribeToLogging
      PubSubMediator.publish(WorkflowExecutionEvent,"workflowevent")
      Thread.sleep(1000L)
      assert(logMessage == "workflowevent")
    }
  }

  "BusinessLogging Actor " must {
    "subscribe a subscriber and pubsubmediator publish null message " in {
      logMessage = "message"
      val actor = TestActorRef(new BusinessLogging() with MockBusinessLogEvent with PubSubMediator)
      actor ! SubscribeToLogging
      PubSubMediator.publish(null,"workflowevent")
      Thread.sleep(1000L)
      assert(logMessage == "message")
    }
  }

}
