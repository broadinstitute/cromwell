package cromwell.core.io

import akka.actor.{Actor, ActorLogging, ActorRef}
import akka.testkit.{TestActorRef, TestProbe}
import common.util.Backoff
import cromwell.core.TestKitSuite
import cromwell.core.io.DefaultIoCommand.DefaultIoSizeCommand
import cromwell.core.path.Path
import cromwell.core.retry.SimpleExponentialBackoff
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.scalatestplus.mockito.MockitoSugar

import scala.concurrent.duration.{FiniteDuration, _}
import scala.language.postfixOps

class IoClientHelperSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with MockitoSugar {

  behavior of "IoClientHelperSpec"

  it should "intercept IoAcks and cancel timers" in {
    val ioActorProbe = TestProbe()
    val delegateProbe = TestProbe()
    val backoff = SimpleExponentialBackoff(100 seconds, 10.hours, 2D, 0D)
    val noResponseTimeout = 3 seconds
    
    val testActor = TestActorRef(new IoClientHelperTestActor(ioActorProbe.ref, delegateProbe.ref, backoff, noResponseTimeout)) 
    
    val command = DefaultIoSizeCommand(mock[Path])
    val response = IoSuccess(command, 5)

    // Send the command
    testActor.underlyingActor.sendMessage(command)

    // Io actor receives the command
    ioActorProbe.expectMsg(command)
    
    // Io actor replies
    ioActorProbe.reply(response)
    
    // delegate should receive the response
    delegateProbe.expectMsg(response)
    
    // And nothing else, meaning the timeout timer has been cancelled
    delegateProbe.expectNoMessage()

    // timeouts map should be empty
    testActor.underlyingActor.timeouts.isEmpty shouldBe true
  }

  it should "intercept IoAcks and cancel timers for a command with context" in {
    val ioActorProbe = TestProbe()
    val delegateProbe = TestProbe()
    val backoff = SimpleExponentialBackoff(100 seconds, 10.hours, 2D, 0D)
    val noResponseTimeout = 3 seconds

    val testActor = TestActorRef(new IoClientHelperTestActor(ioActorProbe.ref, delegateProbe.ref, backoff, noResponseTimeout))

    val commandContext = "context"
    val command = DefaultIoSizeCommand(mock[Path])
    val response = IoSuccess(command, 5)

    // Send the command
    testActor.underlyingActor.sendMessageWithContext(commandContext, command)

    // Io actor receives the command
    ioActorProbe.expectMsg(commandContext -> command)

    // Io actor replies
    ioActorProbe.reply(commandContext -> response)

    // delegate should receive the response
    delegateProbe.expectMsgPF(1 second) {
      case (contextReceived,  responseReceived) if contextReceived == "context" && responseReceived == response =>
    }

    // And nothing else, meaning the timeout timer has been cancelled
    delegateProbe.expectNoMessage()

    // timeouts map should be empty
    testActor.underlyingActor.timeouts.isEmpty shouldBe true
  }

  private case object ServiceUnreachable
  
  private class IoClientHelperTestActor(override val ioActor: ActorRef,
                                delegateTo: ActorRef,
                                backoff: Backoff,
                                noResponseTimeout: FiniteDuration) extends Actor with ActorLogging with IoClientHelper {

    implicit val ioCommandBuilder = DefaultIoCommandBuilder

    override protected def initialBackoff = backoff
    
    context.become(ioReceive orElse receive)

    override def receive: Receive = {
      case message => delegateTo ! message
    }

    def sendMessage(command: IoCommand[_]) = {
      sendIoCommandWithCustomTimeout(command, noResponseTimeout)
    }

    def sendMessageWithContext(context: Any, command: IoCommand[_]) = {
      sendIoCommandWithContext(command, context, noResponseTimeout)
    }

    override protected def onTimeout(message: Any, to: ActorRef): Unit = {
      delegateTo ! ServiceUnreachable
    }
  }

}
