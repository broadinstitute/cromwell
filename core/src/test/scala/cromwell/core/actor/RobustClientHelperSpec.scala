package cromwell.core.actor

import akka.actor.{Actor, ActorLogging, ActorRef}
import akka.testkit.{ImplicitSender, TestActorRef, TestProbe}
import common.util.Backoff
import cromwell.core.TestKitSuite
import cromwell.core.actor.StreamIntegration.BackPressure
import cromwell.core.retry.SimpleExponentialBackoff
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

import scala.concurrent.duration._
import scala.language.postfixOps

class RobustClientHelperSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with ImplicitSender {
  behavior of "RobustClientHelper"
  
  it should "handle Backpressure responses" in {
    val remoteActor = TestProbe()
    val delegateActor = TestProbe()
    
    val margin = 2.second
    val backoff = SimpleExponentialBackoff(1.second, 10.seconds, 2D, 0D)
    val noResponseTimeout = 10 seconds
    val testActor = TestActorRef(new TestActor(delegateActor.ref, backoff, noResponseTimeout))
    
    val messageToSend = TestActor.TestMessage("hello")
    
    //send message
    testActor.underlyingActor.sendMessage(messageToSend, remoteActor.ref)
    
    // remote actor receives message
    remoteActor.expectMsg(messageToSend)
    
    // remote actor sends a backpressure message
    remoteActor.reply(BackPressure(messageToSend))
    
    // remote actor expects request again after backpressureTimeout
    remoteActor.expectMsg(1.second + margin, messageToSend)
    
    // remote actor replies
    remoteActor.reply("world")
    
    // delegate actor receives response
    delegateActor.expectMsg("world")
    
    // remote actor doesn't receives new messages
    remoteActor.expectNoMessage()
    
    // Wait long enough that to make sure that we won't receive a ServiceUnreachable message, meaning the timeout timer
    // has been cancelled. Note that it is the responsibility of the actor to cancel it, the RobustClientHelper does not
    // handle that part.
    delegateActor.expectNoMessage(8 seconds)
  }

  it should "handle a successful response" in {
    val remoteActor = TestProbe()
    val delegateActor = TestProbe()

    val backoff = SimpleExponentialBackoff(1.second, 10.seconds, 2D, 0D)
    val noResponseTimeout = 20 seconds
    val testActor = TestActorRef(new TestActor(delegateActor.ref, backoff, noResponseTimeout))

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
    remoteActor.expectNoMessage()
    delegateActor.expectNoMessage()
  }

  it should "timeout if no response" in {
    val remoteActor = TestProbe()
    val delegateActor = TestProbe()

    val backoff = SimpleExponentialBackoff(1.second, 10.seconds, 2D, 0D)
    val noResponseTimeout = 2 seconds
    val testActor = TestActorRef(new TestActor(delegateActor.ref, backoff, noResponseTimeout))

    val messageToSend = TestActor.TestMessage("hello")

    // send message
    testActor.underlyingActor.sendMessage(messageToSend, remoteActor.ref)

    // remote actor receives message
    remoteActor.expectMsg(messageToSend)

    // remote actor does not reply

    // delegate receives ServiceUnreachable message
    delegateActor.expectMsg(TestActor.ServiceUnreachable)

    // remote actor doesn't receives new messages
    remoteActor.expectNoMessage()
    delegateActor.expectNoMessage()
  }

  it should "reset timeout when backpressured is received" in {
    val remoteActor = TestProbe()
    val delegateActor = TestProbe()
    
    val margin = 1 second
    val backoff = SimpleExponentialBackoff(1.second, 10.seconds, 2D, 0D)
    val noResponseTimeout = 3 seconds
    val testActor = TestActorRef(new TestActor(delegateActor.ref, backoff, noResponseTimeout))

    val messageToSend = TestActor.TestMessage("hello")

    // send message
    testActor.underlyingActor.sendMessage(messageToSend, remoteActor.ref)

    // remote actor receives message
    remoteActor.expectMsg(messageToSend)

    // remote actor sends a backpressure message
    remoteActor.reply(BackPressure(messageToSend))

    // remote actor expects request again after backpressureTimeout
    remoteActor.expectMsg(1.second + margin, messageToSend)

    // remote actor replies
    remoteActor.reply("world")

    // delegate receives ServiceUnreachable message
    delegateActor.expectMsg("world")

    // remote actor doesn't receives new messages
    remoteActor.expectNoMessage()
    // ensure that no ServiceUnreachable message was sent
    delegateActor.expectNoMessage(4 seconds)
  }

  private [actor] object TestActor {
    case class TestMessage(v: String)
    case object ServiceUnreachable
  }
  private class TestActor(delegateTo: ActorRef,
                          backoff: Backoff,
                          noResponseTimeout: FiniteDuration) extends Actor with ActorLogging with RobustClientHelper {

    override def initialBackoff(): Backoff = backoff

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
