package cromwell.services

import akka.actor.PoisonPill
import akka.pattern.ask
import com.typesafe.config.ConfigFactory
import cromwell.core.WorkflowId
import cromwell.services.MetadataServiceActor._

class MetadataServiceSpec extends CromwellServicesSpec {

  val config = ConfigFactory.load("{}")
  val actor = actorSystem.actorOf(MetadataServiceActor.props(config, config))

  val workflowId = WorkflowId.randomId()

  behavior of "MetadataServiceActor"

  /*
  Simple store / retrieve
   */

  val key1 = MetadataKey(workflowId, None, "key1")
  val key2 = MetadataKey(workflowId, None, "key2")

  val event1_1 = MetadataEvent(key1, MetadataValue("value1"))
  val event1_2 = MetadataEvent(key1, MetadataValue("value2"))
  val event2_1 = MetadataEvent(key2, MetadataValue("value1"))

  it should "Store values for different keys" in {
    val putAction1 = PutMetadataAction(event1_1)
    val putAction2 = PutMetadataAction(event1_2)
    val putAction3 = PutMetadataAction(event2_1)
    (for {
      response1 <- (actor ? putAction1).mapTo[MetadataServiceResponse]
      response2 <- (actor ? putAction2).mapTo[MetadataServiceResponse]
      response3 <- (actor ? putAction3).mapTo[MetadataServiceResponse]
      _ = response1 shouldBe MetadataPutAcknowledgement(putAction1)
      _ = response2 shouldBe MetadataPutAcknowledgement(putAction2)
      _ = response3 shouldBe MetadataPutAcknowledgement(putAction3)
    } yield ()).futureValue
  }

  it should "Retrieve the correct values for different keys" in {
    val query1 = MetadataQuery.forKey(key1)
    val query2 = MetadataQuery.forKey(key2)

    (for {
      response1 <- (actor ? GetMetadataQueryAction(query1)).mapTo[MetadataServiceResponse]
      _ = response1 shouldBe MetadataLookupResponse(query1, Seq(event1_1, event1_2))

      response2 <- (actor ? GetMetadataQueryAction(query2)).mapTo[MetadataServiceResponse]
      _ = response2 shouldBe MetadataLookupResponse(query2, Seq(event2_1))
    } yield ()).futureValue
  }

  override def afterAll = {
    actor ! PoisonPill
  }
}
