package cromwell.services.keyvalue

import java.io.IOException

import akka.actor.{Actor, ActorLogging, ActorRef}
import akka.testkit.{TestActorRef, TestProbe}
import cromwell.core.TestKitSuite
import cromwell.services.keyvalue.KeyValueServiceActor._
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

import scala.concurrent.{Await, ExecutionContextExecutor}
import scala.concurrent.duration._
import scala.language.postfixOps

class KvClientSpec extends TestKitSuite with AnyFlatSpecLike with Matchers {

  implicit val ec: ExecutionContextExecutor = system.dispatcher

  behavior of "KvClient"

  it should "Correctly forward multiple requests and responses" in {
    val serviceActorProbe = TestProbe("serviceActorProbe")
    val kvTestClient = TestActorRef(factory = new KvTestClientActor(serviceActorProbe.ref), name = "kvTestClient")

    val scopedKey1 = ScopedKey(null, null, "key1")
    val scopedKey2 = ScopedKey(null, null, "key2")
    val putRequest = KvPut(KvPair(scopedKey1, "value1"))
    val getRequest = KvGet(scopedKey2)
    val putResponse = KvFailure(putRequest, new IOException())
    val getResponse = KvPair(scopedKey2, "value2")

    val requests = Seq(putRequest, getRequest)
    val futureResult = kvTestClient.underlyingActor.makeKvRequest(requests)

    serviceActorProbe.expectMsgAllOf(putRequest, getRequest)
    serviceActorProbe.expectNoMessage(max = 50 milliseconds)

    kvTestClient.underlyingActor.currentKvClientRequests.size should be(2)

    kvTestClient.tell(getResponse, sender = serviceActorProbe.ref)
    serviceActorProbe.expectNoMessage(max = 50 milliseconds)
    futureResult.isCompleted should be(false)

    kvTestClient.underlyingActor.currentKvClientRequests.size should be(1)

    kvTestClient.tell(putResponse, sender = serviceActorProbe.ref)
    serviceActorProbe.expectNoMessage(max = 50 milliseconds)

    // Make sure the future completes promptly and the original order is preserved:
    Await.result(futureResult, atMost = 100 milliseconds) should be(Seq(putResponse, getResponse))
    serviceActorProbe.expectNoMessage(max = 50 milliseconds)

    kvTestClient.underlyingActor.currentKvClientRequests.size should be(0)
  }
}

class KvTestClientActor(val serviceRegistryActor: ActorRef) extends Actor with ActorLogging with KvClient {
  override def receive: Receive = kvClientReceive orElse Actor.ignoringBehavior
}


