package cromwell.services.metadata

import java.util.UUID

import akka.actor.{ActorRef, Props}
import akka.testkit.TestProbe
import cromwell.core.{TestKitSuite, WorkflowId}
import cromwell.services.metadata.MetadataService.{BuildMetadataJsonAction, GetMetadataAction}
import cromwell.services.metadata.impl.ReadMetadataRegulatorActor
import cromwell.services.{MetadataJsonResponse, SuccessfulMetadataJsonResponse}
import org.scalatest.flatspec.AnyFlatSpecLike
import spray.json.JsObject

import scala.concurrent.duration._

class ReadMetadataRegulatorActorSpec extends TestKitSuite with AnyFlatSpecLike {

  behavior of "ReadMetadataRegulatorActor"

  val defaultTimeout: FiniteDuration = 1.second

  def readMetadataRegulatorActorTestProps(metadataBuilderActors: Iterator[ActorRef] = Iterator.empty,
                                          readMetadataWorkers: Iterator[ActorRef] = Iterator.empty): Props = {
    val metadataBuilderActorProps = () => ???
    val readMetadataWorkerProps = () => ???
    Props {
      new ReadMetadataRegulatorActor(metadataBuilderActorProps, readMetadataWorkerProps) {
        override def createMetadataBuilderActor(workflowId: WorkflowId): ActorRef = metadataBuilderActors.next()

        override def createReadMetadataWorker(): ActorRef = readMetadataWorkers.next()
      }
    }
  }

  def createRequest: BuildMetadataJsonAction =
    GetMetadataAction(MetadataQuery(WorkflowId.randomId(), None, None, None, None, expandSubWorkflows = false))

  def createResponse(action: BuildMetadataJsonAction): MetadataJsonResponse =
    SuccessfulMetadataJsonResponse(action, JsObject.empty)

  it should "respond to a single request" in {

    val sender = TestProbe("actionSender")
    val mockMetadataBuilderActor = TestProbe("mockMetadataBuilderActor")

    val readMetadataRegulatorActor =
      system.actorOf(
        props = readMetadataRegulatorActorTestProps(Iterator.single(mockMetadataBuilderActor.ref)),
        name = s"ReadMetadataRegulatorActor-${UUID.randomUUID()}")

    val request = createRequest
    val response = createResponse(request)

    sender.send(readMetadataRegulatorActor, request)

    mockMetadataBuilderActor.expectMsg(defaultTimeout, request)
    mockMetadataBuilderActor.reply(response)

    sender.expectMsg(defaultTimeout, response)
  }

  it should "respond to different requests from same sender" in {

    val sender = TestProbe("actionSender")
    val mockMetadataBuilderActor1 = TestProbe("mockMetadataBuilderActor1")
    val mockMetadataBuilderActor2 = TestProbe("mockMetadataBuilderActor2")

    val readMetadataRegulatorActor =
      system.actorOf(
        props = readMetadataRegulatorActorTestProps(Iterator(mockMetadataBuilderActor1.ref, mockMetadataBuilderActor2.ref)),
        name = s"ReadMetadataRegulatorActor-${UUID.randomUUID()}")

    val request1 = createRequest
    val request2 = createRequest

    val response1 = createResponse(request1)
    val response2 = createResponse(request2)

    sender.send(readMetadataRegulatorActor, request1)
    sender.send(readMetadataRegulatorActor, request2)

    mockMetadataBuilderActor1.expectMsg(defaultTimeout, request1)
    mockMetadataBuilderActor1.reply(response1)

    mockMetadataBuilderActor2.expectMsg(defaultTimeout, request2)
    mockMetadataBuilderActor2.reply(response2)

    sender.expectMsg(defaultTimeout, response1)
    sender.expectMsg(defaultTimeout, response2)
  }

  it should "respond to different requests from different senders" in {

    val sender1 = TestProbe("actionSender1")
    val sender2 = TestProbe("actionSender2")

    val mockMetadataBuilderActor1 = TestProbe("mockMetadataBuilderActor1")
    val mockMetadataBuilderActor2 = TestProbe("mockMetadataBuilderActor2")

    val readMetadataRegulatorActor =
      system.actorOf(
        props = readMetadataRegulatorActorTestProps(Iterator(mockMetadataBuilderActor1.ref, mockMetadataBuilderActor2.ref)),
        name = s"ReadMetadataRegulatorActor-${UUID.randomUUID()}")

    val request1 = createRequest
    val request2 = createRequest

    val response1 = createResponse(request1)
    val response2 = createResponse(request2)

    sender1.send(readMetadataRegulatorActor, request1)
    sender2.send(readMetadataRegulatorActor, request2)

    mockMetadataBuilderActor1.expectMsg(defaultTimeout, request1)
    mockMetadataBuilderActor1.reply(response1)

    mockMetadataBuilderActor2.expectMsg(defaultTimeout, request2)
    mockMetadataBuilderActor2.reply(response2)

    sender1.expectMsg(defaultTimeout, response1)
    sender2.expectMsg(defaultTimeout, response2)
  }

  it should "respond to same request from different senders" in {

    val sender1 = TestProbe("actionSender1")
    val sender2 = TestProbe("actionSender2")

    val mockMetadataBuilderActor = TestProbe("mockMetadataBuilderActor")

    val readMetadataRegulatorActor =
      system.actorOf(
        props = readMetadataRegulatorActorTestProps(Iterator.single(mockMetadataBuilderActor.ref)),
        name = s"ReadMetadataRegulatorActor-${UUID.randomUUID()}")

    val request = createRequest
    val response = createResponse(request)

    sender1.send(readMetadataRegulatorActor, request)
    sender2.send(readMetadataRegulatorActor, request)

    mockMetadataBuilderActor.expectMsg(defaultTimeout, request)
    mockMetadataBuilderActor.reply(response)

    sender1.expectMsg(defaultTimeout, response)
    sender2.expectMsg(defaultTimeout, response)
  }

  it should "respond to same requests from same sender" in {

    val sender = TestProbe("actionSender")
    val mockMetadataBuilderActor = TestProbe("mockMetadataBuilderActor")

    val readMetadataRegulatorActor =
      system.actorOf(
        props = readMetadataRegulatorActorTestProps(Iterator.single(mockMetadataBuilderActor.ref)),
        name = s"ReadMetadataRegulatorActor-${UUID.randomUUID()}")

    val request = createRequest
    val response = createResponse(request)

    sender.send(readMetadataRegulatorActor, request)
    sender.send(readMetadataRegulatorActor, request)

    mockMetadataBuilderActor.expectMsg(defaultTimeout, request)
    mockMetadataBuilderActor.reply(response)

    sender.expectMsg(defaultTimeout, response)
    sender.expectMsg(defaultTimeout, response)
  }

}
