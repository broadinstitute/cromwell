package cromwell.core.actor

import akka.actor.{Actor, ActorLogging, ActorRef}
import akka.testkit.{ImplicitSender, TestActorRef, TestProbe}
import cromwell.core.TestKitSuite
import cromwell.core.actor.StreamIntegration.Backpressure
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration._
import scala.language.postfixOps

class RobustClientHelperSpec extends TestKitSuite with FlatSpecLike with Matchers with ImplicitSender {
  behavior of "RobustClientHelper"
  
  it should "handle Backpressure responses" in {
    val remoteActor = TestProbe()
    val delegateActor = TestProbe()
    
    val margin = 1 second
    val backpressureTimeout = 1 second
    val noResponseTimeout = 10 seconds
    val testActor = TestActorRef(new TestActor(delegateActor.ref, backpressureTimeout, noResponseTimeout))
    
    val messageToSend = TestActor.TestMessage("hello")
    
    //send message
    testActor.underlyingActor.sendMessage(messageToSend, remoteActor.ref)
    
    // remote actor receives message
    remoteActor.expectMsg(messageToSend)
    
    // remote actor sends a backpressure message
    remoteActor.reply(Backpressure(messageToSend))
    
    // remote actor expects request again after backpressureTimeout
    remoteActor.expectMsg(backpressureTimeout + margin, messageToSend)
    
    // remote actor replies
    remoteActor.reply("world")
    
    // delegate actor receives response
    delegateActor.expectMsg("world")
    
    // remote actor doesn't receives new messages
    remoteActor.expectNoMsg()
    
    // Wait long enough that to make sure that we won't receive a ServiceUnreachable message, meaning the timeout timer
    // has been cancelled. Note that it is the responsibility of the actor to cancel it, the RobustClientHelper does not
    // handle that part.
    delegateActor.expectNoMsg(8 seconds)
  }

  it should "handle a successful response" in {
    val remoteActor = TestProbe()
    val delegateActor = TestProbe()

    val backpressureTimeout = 1 second
    val noResponseTimeout = 20 seconds
    val testActor = TestActorRef(new TestActor(delegateActor.ref, backpressureTimeout, noResponseTimeout))

    val messageToSend = TestActor.TestMessage("hello")

    // send message
    testActor.underlyingActor.sendMessage(messageToSend, remoteActor.ref)

    // remote actor receives message
    remoteActor.expectMsg(messageToSend)
    
    // remote actor replies
    remoteActor.reply("world")
    
    // delegate receives response
    delegateActor.expectMsg("world")
    
    // remote actor doesn't receives new messages
    remoteActor.expectNoMsg()
    delegateActor.expectNoMsg()
  }

  it should "timeout if no response" in {
    val remoteActor = TestProbe()
    val delegateActor = TestProbe()

    val backpressureTimeout = 1 second
    val noResponseTimeout = 2 seconds
    val testActor = TestActorRef(new TestActor(delegateActor.ref, backpressureTimeout, noResponseTimeout))

    val messageToSend = TestActor.TestMessage("hello")

    // send message
    testActor.underlyingActor.sendMessage(messageToSend, remoteActor.ref)

    // remote actor receives message
    remoteActor.expectMsg(messageToSend)

    // remote actor does not reply

    // delegate receives ServiceUnreachable message
    delegateActor.expectMsg(TestActor.ServiceUnreachable)

    // remote actor doesn't receives new messages
    remoteActor.expectNoMsg()
    delegateActor.expectNoMsg()
  }

  it should "reset timeout when backpressured is received" in {
    val remoteActor = TestProbe()
    val delegateActor = TestProbe()
    
    val margin = 1 second
    val backpressureTimeout = 1 second
    val noResponseTimeout = 3 seconds
    val testActor = TestActorRef(new TestActor(delegateActor.ref, backpressureTimeout, noResponseTimeout))

    val messageToSend = TestActor.TestMessage("hello")

    // send message
    testActor.underlyingActor.sendMessage(messageToSend, remoteActor.ref)

    // remote actor receives message
    remoteActor.expectMsg(messageToSend)

    // remote actor sends a backpressure message
    remoteActor.reply(Backpressure(messageToSend))

    // remote actor expects request again after backpressureTimeout
    remoteActor.expectMsg(backpressureTimeout + margin, messageToSend)

    // remote actor replies
    remoteActor.reply("world")

    // delegate receives ServiceUnreachable message
    delegateActor.expectMsg("world")

    // remote actor doesn't receives new messages
    remoteActor.expectNoMsg()
    // ensure that no ServiceUnreachable message was sent
    delegateActor.expectNoMsg(4 seconds)
  }
  
  it should "randomize backpressure timings" in {
    val delegateActor = TestProbe()
    val backpressureTimeout = 100 seconds
    val noResponseTimeout = 3 seconds
    val randomizeFactor = 0.5D
    
    val testActor = TestActorRef(new TestActor(delegateActor.ref, backpressureTimeout, noResponseTimeout, randomizeFactor)).underlyingActor
    
    val randomBackpressures = 0 until 10 map { _ =>
      val time = testActor.generateBackpressureTime
      time.gt(50.seconds) shouldBe true
      time.lt(150.seconds) shouldBe true
      time
    }
    
    // They should all be different
    randomBackpressures.distinct.size > 1 shouldBe true
  }
  
  private [actor] object TestActor {
    case class TestMessage(v: String)
    case object ServiceUnreachable
  }
  private class TestActor(delegateTo: ActorRef,
                          override val backpressureTimeout: FiniteDuration,
                          noResponseTimeout: FiniteDuration,
                          override val backpressureRandomizerFactor: Double = 0.5D) extends Actor with ActorLogging with RobustClientHelper {
    
    context.become(robustReceive orElse receive)
    var messageSent: Any = _
    
    override def receive: Receive = {
      case message =>
        cancelTimeout(messageSent)
        delegateTo ! message
    }
    
    def sendMessage(message: Any, to: ActorRef) = {
      messageSent = message
      robustSend(message, to, noResponseTimeout)
    }

    override protected def onTimeout(message: Any, to: ActorRef): Unit = {
      delegateTo ! TestActor.ServiceUnreachable
    }
  }
}
