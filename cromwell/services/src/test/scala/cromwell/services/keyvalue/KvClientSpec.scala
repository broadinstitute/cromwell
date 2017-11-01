package cromwell.services.keyvalue

import java.io.IOException

import akka.actor.{Actor, ActorLogging, ActorRef, ActorSystem}
import akka.testkit.{TestActorRef, TestKit, TestProbe}
import cromwell.services.keyvalue.KeyValueServiceActor._
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.Await
import scala.concurrent.duration._
import scala.language.postfixOps

class KvClientSpec extends TestKit(ActorSystem("KvClientSpec")) with FlatSpecLike with Matchers {

  implicit val ec = system.dispatcher

  behavior of "KvClient"

  it should "Correctly forward multiple requests and responses" in {
    val serviceActorProbe = TestProbe()
    val kvTestClient = TestActorRef(new KvTestClientActor(serviceActorProbe.ref))

    val scopedKey1 = ScopedKey(null, null, "key1")
    val scopedKey2 = ScopedKey(null, null, "key2")
    val putRequest = KvPut(KvPair(scopedKey1, Some("value1")))
    val getRequest = KvGet(scopedKey2)
    val putResponse = KvFailure(putRequest, new IOException())
    val getResponse = KvPair(scopedKey2, Some("value2"))

    val requests = Seq(putRequest, getRequest)
    val futureResult = kvTestClient.underlyingActor.makeKvRequest(requests)

    serviceActorProbe.expectMsgAllOf(putRequest, getRequest)
    serviceActorProbe.expectNoMsg(max = 50 milliseconds)

    kvTestClient.underlyingActor.currentKvClientRequests.size should be(2)

    kvTestClient.tell(getResponse, sender = serviceActorProbe.ref)
    serviceActorProbe.expectNoMsg(max = 50 milliseconds)
    futureResult.isCompleted should be(false)

    kvTestClient.underlyingActor.currentKvClientRequests.size should be(1)

    kvTestClient.tell(putResponse, sender = serviceActorProbe.ref)
    serviceActorProbe.expectNoMsg(max = 50 milliseconds)

    // Make sure the future completes promptly and the original order is preserved:
    Await.result(futureResult, atMost = 100 milliseconds) should be(Seq(putResponse, getResponse))
    serviceActorProbe.expectNoMsg(max = 50 milliseconds)

    kvTestClient.underlyingActor.currentKvClientRequests.size should be(0)
  }
}

class KvTestClientActor(val serviceRegistryActor: ActorRef) extends Actor with ActorLogging with KvClient {
  override def receive = kvClientReceive orElse Actor.ignoringBehavior
}


