package cromwell.services.metadata.impl

import java.time.OffsetDateTime

import akka.pattern._
import com.typesafe.config.ConfigFactory
import cromwell.core.WorkflowId
import cromwell.services.ServicesSpec
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata._
import org.scalatest.concurrent.Eventually._
import org.scalatest.concurrent.PatienceConfiguration.{Interval, Timeout}

import scala.concurrent.duration._

class MetadataServiceActorSpec extends ServicesSpec("Metadata") {
  import MetadataServiceActorSpec.Config
  val config = ConfigFactory.parseString(Config)
  val actor = system.actorOf(MetadataServiceActor.props(config, config))

    val workflowId = WorkflowId.randomId()

    /*
    Simple store / retrieve
     */

    val key1 = MetadataKey(workflowId, None, "key1")
    val key2 = MetadataKey(workflowId, None, "key2")
    val supJob = MetadataJobKey("sup.sup", None, 1)
    val key3 = MetadataKey(workflowId, Option(supJob), "dog")
    val moment = OffsetDateTime.now

    val event1_1 = MetadataEvent(key1, Option(MetadataValue("value1")), moment)
    val event1_2 = MetadataEvent(key1, Option(MetadataValue("value2")), moment)
    val event2_1 = MetadataEvent(key2, Option(MetadataValue("value1")), moment)
    val event3_1 = MetadataEvent(key3, Option(MetadataValue("value3")), moment)
    val event3_2 = MetadataEvent(key3, None, moment)

  "MetadataServiceActor" should {
    "Store values for different keys and then retrieve those values" in {
      val putAction1 = PutMetadataAction(event1_1)
      val putAction2 = PutMetadataAction(event1_2)
      val putAction3 = PutMetadataAction(event2_1, event3_1, event3_2)
      actor ! putAction1
      actor ! putAction2
      actor ! putAction3

      val query1 = MetadataQuery.forKey(key1)
      val query2 = MetadataQuery.forKey(key2)
      val query3 = MetadataQuery.forKey(key3)
      val query4 = MetadataQuery.forWorkflow(workflowId)
      val query5 = MetadataQuery.forJob(workflowId, supJob)

      eventually(Timeout(10.seconds), Interval(2.seconds)) {
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
    }
  }
}

object MetadataServiceActorSpec {
  val Config =
    """
      |services.MetadataService.db-batch-size = 3
      |services.MetadataService.db-flush-rate = 100 millis
    """.stripMargin
}
