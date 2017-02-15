package cromwell.core.actor

import akka.actor.{Actor, ActorLogging, ActorRef, Props}
import akka.stream.QueueOfferResult.Dropped
import akka.stream.scaladsl.Source
import akka.stream.{ActorMaterializer, OverflowStrategy}
import akka.testkit.{ImplicitSender, TestActorRef}
import cromwell.core.TestKitSuite
import cromwell.core.actor.StreamIntegration._
import cromwell.core.actor.TestStreamActor.{TestStreamActorCommand, TestStreamActorContext}
import org.scalatest.{FlatSpecLike, Matchers}

class StreamActorHelperSpec extends TestKitSuite with FlatSpecLike with Matchers with ImplicitSender {
  behavior of "StreamActorHelper"
  
  implicit val materializer = ActorMaterializer()

  it should "catch EnqueueResponse message" in {
    val actor = TestActorRef(Props(new TestStreamActor(1)))
    val command = new TestStreamActorCommand
    actor ! command
    expectMsg("hello")
    system stop actor
  }
  
  it should "send a backpressure message when messages are dropped by the queue" in {
    val actor = TestActorRef(new TestStreamActor(1))
    val command = new TestStreamActorCommand

    actor ! EnqueueResponse(Dropped, TestStreamActorContext(command, self))

    expectMsg(Backpressure(command))

    system stop actor
  }
}


private object TestStreamActor {
  class TestStreamActorCommand
  case class TestStreamActorContext(request: TestStreamActorCommand, replyTo: ActorRef) extends StreamContext
}

private class TestStreamActor(queueSize: Int)(implicit override val materializer: ActorMaterializer) extends Actor with ActorLogging with StreamActorHelper[TestStreamActorContext] {
  
  override protected def actorReceive: Receive = {
    case command: TestStreamActorCommand =>
      val replyTo = sender()
      val commandContext = TestStreamActorContext(command, replyTo)
      sendToStream(commandContext)
  }

  override protected val streamSource = Source.queue[TestStreamActorContext](queueSize, OverflowStrategy.dropNew)
    .map{ ("hello", _) }
}
