package cromwell.services

import akka.actor.{Actor, ActorRef, Props}
import akka.testkit.{ImplicitSender, TestActorRef, TestProbe}
import common.assertion.ManyTimes.intWithTimes
import cromwell.core.TestKitSuite
import cromwell.services.IoActorRequesterSpec._
import cromwell.services.ServiceRegistryActor.{IoActorRef, NoIoActorRefAvailable, RequestIoActorRef}
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

import scala.concurrent.Future
import scala.concurrent.duration._
import scala.util.{Failure, Success}

class IoActorRequesterSpec extends TestKitSuite with ImplicitSender with AnyFlatSpecLike with Matchers with Eventually {

  behavior of "IoActorRequester"

  it should "Allow actors to request and receive a reference to the IO Actor" in {
    val serviceRegistry = TestProbe()
    val ioActor = TestProbe()
    val requester = TestActorRef[SimpleIoActorRequester](simpleIoActorRequesterProps(serviceRegistry.ref), "simple-io-requester-request-receive")

    serviceRegistry.expectMsg(RequestIoActorRef)
    serviceRegistry.lastSender ! IoActorRef(ioActor.ref)

    eventually {
      requester.underlyingActor.ioActor.value should be(Some(Success(ioActor.ref)))
    }

    requester.stop()
  }

  it should "keep asking if the IoActor is not available immediately" in {
    val serviceRegistry = TestProbe()
    val ioActor = TestProbe()
    val requester = TestActorRef[SimpleIoActorRequester](simpleIoActorRequesterProps(serviceRegistry.ref), "simple-io-requester-not-ready-resilience")

    10.times {
      serviceRegistry.expectMsg(RequestIoActorRef)
      serviceRegistry.lastSender ! NoIoActorRefAvailable
    }

    serviceRegistry.expectMsg(RequestIoActorRef)
    serviceRegistry.lastSender ! IoActorRef(ioActor.ref)

    eventually {
      requester.underlyingActor.ioActor.value should be(Some(Success(ioActor.ref)))
    }

    requester.stop()
  }

  it should "fail if the IoActorRequest receives an unexpected response" in {
    val serviceRegistry = TestProbe()
    val requester = TestActorRef[SimpleIoActorRequester](simpleIoActorRequesterProps(serviceRegistry.ref), "simple-io-requester-not-ready-resilience")

    serviceRegistry.expectMsg(RequestIoActorRef)
    serviceRegistry.lastSender ! "!!!Bad Response!!!"

    eventually {
      requester.underlyingActor.ioActor.value match {
        case Some(Failure(exception)) => exception.getMessage should startWith("Programmer Error: Unexpected response to a RequestIoActor message")
        case Some(Success(other)) => fail(s"Expected the IoActor ref to be empty, but got: $other")
        case None => fail("No IoActor received")
      }
    }

    requester.stop()
  }
}

object IoActorRequesterSpec {
  class SimpleIoActorRequester(override val serviceRegistryActor: ActorRef) extends Actor with IoActorRequester {

    val ioActor: Future[ActorRef] = requestIoActor(backoffInterval = 10.millis)

    override def receive: Receive = {
      case other => throw new RuntimeException(s"Didn't expect a message but got: $other")
    }
  }

  def simpleIoActorRequesterProps(serviceRegistryActor: ActorRef): Props = Props(new SimpleIoActorRequester(serviceRegistryActor))
}
