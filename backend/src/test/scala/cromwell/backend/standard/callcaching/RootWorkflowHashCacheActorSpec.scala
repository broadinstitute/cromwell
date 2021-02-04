package cromwell.backend.standard.callcaching

import akka.actor.Props
import akka.testkit._
import cromwell.backend.standard.callcaching.RootWorkflowFileHashCacheActor.{IoHashCommandWithContext, _}
import cromwell.core.actor.RobustClientHelper.RequestTimeout
import cromwell.core.callcaching.HashKey
import cromwell.core.io.DefaultIoCommand.DefaultIoHashCommand
import cromwell.core.io.IoSuccess
import cromwell.core.path.DefaultPathBuilder
import cromwell.core.{TestKitSuite, WorkflowId}
import org.scalatest.flatspec.AnyFlatSpecLike

import scala.concurrent.duration._

class RootWorkflowHashCacheActorSpec extends TestKitSuite with ImplicitSender
  with AnyFlatSpecLike {

  private val fakeWorkflowId = WorkflowId.randomId()
  private val fakeFileName = "fakeFileName"

  it should "properly handle the situation when response from IoActor received after the timeout" in {
    val ioActorProbe = TestProbe("ioActorProbe-without-timer")
    val rootWorkflowFileHashCacheActor = system.actorOf(
      props = Props(new RootWorkflowFileHashCacheActor(ioActorProbe.ref, fakeWorkflowId) {
        override lazy val defaultIoTimeout: FiniteDuration = 1.second
      }),
      name = "rootWorkflowFileHashCacheActor-without-timer",
    )

    val ioHashCommandWithContext = IoHashCommandWithContext(DefaultIoHashCommand(DefaultPathBuilder.build("").get), FileHashContext(HashKey(checkForHitOrMiss = false, List.empty), fakeFileName))
    rootWorkflowFileHashCacheActor ! ioHashCommandWithContext

    //wait for timeout
    Thread.sleep(2000)

    EventFilter.info(msgIoAckWithNoRequesters.format(fakeFileName), occurrences = 1).intercept {
      ioActorProbe.send(rootWorkflowFileHashCacheActor, ioHashCommandWithContext.fileHashContext -> IoSuccess(ioHashCommandWithContext.ioHashCommand, "Successful result"))
    }
  }

  it should "properly handle the situation when timeout occurs when response from IoActor has already been received, but timer has not yet been disabled" in {
    val ioActorProbe = TestProbe("ioActorProbe-with-timer")
    val rootWorkflowFileHashCacheActor = system.actorOf(
      Props(new RootWorkflowFileHashCacheActor(ioActorProbe.ref, fakeWorkflowId) {
        // Effectively disabling automatic timeout firing here. We'll send RequestTimeout ourselves
        override lazy val defaultIoTimeout: FiniteDuration = 1.hour
      }),
      "rootWorkflowFileHashCacheActor-with-timer",
    )

    val ioHashCommandWithContext = IoHashCommandWithContext(DefaultIoHashCommand(DefaultPathBuilder.build("").get), FileHashContext(HashKey(checkForHitOrMiss = false, List.empty), fakeFileName))
    rootWorkflowFileHashCacheActor ! ioHashCommandWithContext

    val hashVal = "Success"
    EventFilter.info(msgTimeoutAfterIoAck.format(s"FileHashSuccess($hashVal)", ioHashCommandWithContext.fileHashContext.file), occurrences = 1).intercept {
      ioActorProbe.send(rootWorkflowFileHashCacheActor, (ioHashCommandWithContext.fileHashContext, IoSuccess(ioHashCommandWithContext.ioHashCommand, hashVal)))
      Thread.sleep(2000) // wait for actor to put value into cache
      ioActorProbe.send(rootWorkflowFileHashCacheActor, RequestTimeout(ioHashCommandWithContext.fileHashContext -> ioHashCommandWithContext.ioHashCommand, rootWorkflowFileHashCacheActor))
    }
  }
}
