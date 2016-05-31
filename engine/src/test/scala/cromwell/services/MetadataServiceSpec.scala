package cromwell.services

import java.time.OffsetDateTime

import akka.actor.PoisonPill
import akka.pattern.ask
import com.typesafe.config.ConfigFactory
import cromwell.core.WorkflowId
import cromwell.services.MetadataServiceActor._

class MetadataServiceSpec extends CromwellServicesSpec {

  val config = ConfigFactory.load("{}")
  // TODO: PBE: No EngineMetadataServiceActor.props until circular dependencies fixed.
  //val actor = actorSystem.actorOf(MetadataServiceActor.props(config, config))
  val actor = actorSystem.actorOf(akka.actor.Props(new EngineMetadataServiceActor(config, config)))

  val workflowId = WorkflowId.randomId()

  behavior of "MetadataServiceActor"

  /*
  Simple store / retrieve
   */

  val key1 = MetadataKey(workflowId, None, "key1")
  val key2 = MetadataKey(workflowId, None, "key2")
  val supJob = MetadataJobKey("sup.sup", None, 1)
  val key3 = MetadataKey(workflowId, Option(supJob), "dog")
  val moment = OffsetDateTime.now

  val event1_1 = MetadataEvent(key1, MetadataValue("value1"), moment)
  val event1_2 = MetadataEvent(key1, MetadataValue("value2"), moment)
  val event2_1 = MetadataEvent(key2, MetadataValue("value1"), moment)
  val event3_1 = MetadataEvent(key3, MetadataValue("value3"), moment)
  val event3_2 = MetadataEvent(key3, MetadataValue(null), moment)

  it should "Store values for different keys" in {
    val putAction1 = PutMetadataAction(event1_1)
    val putAction2 = PutMetadataAction(event1_2)
    val putAction3 = PutMetadataAction(event2_1)
    val putAction4 = PutMetadataAction(event3_1)
    val putAction5 = PutMetadataAction(event3_2)
    (for {
      response1 <- (actor ? putAction1).mapTo[MetadataServiceResponse]
      response2 <- (actor ? putAction2).mapTo[MetadataServiceResponse]
      response3 <- (actor ? putAction3).mapTo[MetadataServiceResponse]
      response4 <- (actor ? putAction4).mapTo[MetadataServiceResponse]
      response5 <- (actor ? putAction5).mapTo[MetadataServiceResponse]
      _ = response1 shouldBe MetadataPutAcknowledgement(putAction1)
      _ = response2 shouldBe MetadataPutAcknowledgement(putAction2)
      _ = response3 shouldBe MetadataPutAcknowledgement(putAction3)
      _ = response4 shouldBe MetadataPutAcknowledgement(putAction4)
      _ = response5 shouldBe MetadataPutAcknowledgement(putAction5)
    } yield ()).futureValue
  }

  it should "Retrieve the correct values for different keys" in {
    val query1 = MetadataQuery.forKey(key1)
    val query2 = MetadataQuery.forKey(key2)
    val query3 = MetadataQuery.forKey(key3)
    val query4 = MetadataQuery.forWorkflow(workflowId)
    val query5 = MetadataQuery.forJob(workflowId, supJob)

    (for {
      response1 <- (actor ? GetMetadataQueryAction(query1)).mapTo[MetadataServiceResponse]
      _ = response1 shouldBe MetadataLookupResponse(query1, Seq(event1_1, event1_2))

      response2 <- (actor ? GetMetadataQueryAction(query2)).mapTo[MetadataServiceResponse]
      _ = response2 shouldBe MetadataLookupResponse(query2, Seq(event2_1))

      response3 <- (actor ? GetMetadataQueryAction(query3)).mapTo[MetadataServiceResponse]
      _ = response3 shouldBe MetadataLookupResponse(query3, Seq(event3_1, event3_2))

      response4 <- (actor ? GetMetadataQueryAction(query4)).mapTo[MetadataServiceResponse]
      _ = response4 shouldBe MetadataLookupResponse(query4, Seq(event1_1, event1_2, event2_1, event3_1, event3_2))

      response5 <- (actor ? GetMetadataQueryAction(query5)).mapTo[MetadataServiceResponse]
      _ = response5 shouldBe MetadataLookupResponse(query5, Seq(event3_1, event3_2))

    } yield ()).futureValue
  }

  override def afterAll() = {
    actor ! PoisonPill
  }
}
