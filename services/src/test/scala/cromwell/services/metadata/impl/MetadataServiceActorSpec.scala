package cromwell.services.metadata.impl

import java.time.OffsetDateTime

import akka.pattern._
import akka.testkit.TestProbe
import com.typesafe.config.ConfigFactory
import cromwell.core._
import cromwell.services.ServicesSpec
import cromwell.services.metadata.MetadataService._
import cromwell.services.metadata._
import cromwell.services.metadata.impl.builder.MetadataBuilderActor.BuiltMetadataResponse

import scala.concurrent.Await
import org.scalatest.concurrent.Eventually._
import org.scalatest.concurrent.PatienceConfiguration.{Interval, Timeout}

import spray.json._

import scala.concurrent.duration._

class MetadataServiceActorSpec extends ServicesSpec("Metadata") {
  import MetadataServiceActorSpec.Config
  val config = ConfigFactory.parseString(Config)
  val actor = system.actorOf(MetadataServiceActor.props(config, config, TestProbe().ref), "MetadataServiceActor-for-MetadataServiceActorSpec")

    val workflowId = WorkflowId.randomId()

    /*
    Simple store / retrieve
     */

    val key1 = MetadataKey(workflowId, None, "key1")
    val key2 = MetadataKey(workflowId, None, "key2")
    val supJob = MetadataJobKey("sup.sup", None, 1)
    val key3 = MetadataKey(workflowId, Option(supJob), "dog")
  val moment = OffsetDateTime.now.minusMinutes(1)

  val event1_1 = MetadataEvent(key1, Option(MetadataValue("value1")), moment.plusSeconds(1))
  val event1_2 = MetadataEvent(key1, Option(MetadataValue("value2")), moment.plusSeconds(2))
  val event2_1 = MetadataEvent(key2, Option(MetadataValue("value1")), moment.plusSeconds(3))
  val event3_1 = MetadataEvent(key3, Option(MetadataValue("value3")), moment.plusSeconds(4))
  val event3_2 = MetadataEvent(key3, None, moment.plusSeconds(5))

  override def beforeAll: Unit = {

    // Even though event1_1 arrives second, the older timestamp should mean it does not replace event1_2:
    val putAction2 = PutMetadataAction(event1_2)
    val putAction1 = PutMetadataAction(event1_1)
    val putAction3 = PutMetadataAction(event2_1, event3_1, event3_2)

    actor ! putAction1
    actor ! putAction2
    actor ! putAction3
  }

  val query1 = MetadataQuery.forKey(key1)
  val query2 = MetadataQuery.forKey(key2)
  val query3 = MetadataQuery.forKey(key3)
  val query4 = MetadataQuery.forWorkflow(workflowId)
  val query5 = MetadataQuery.forJob(workflowId, supJob)

  def expectConstructedMetadata(query: MetadataQuery, expectation: String) = {

  }

  val testCases = List[(String, MetadataQuery, String)] (
    ("query1", query1, s"""{
                          |  "key1": "value2",
                          |  "calls": {},
                          |  "id": "$workflowId"
                          |}""".stripMargin),
    ("query2", query2, s"""{
                          |  "key2": "value1",
                          |  "calls": {},
                          |  "id": "$workflowId"
                          |}""".stripMargin),
    ("query3", query3, s"""{
                          |  "calls": {
                          |    "sup.sup": [{
                          |      "dog": {},
                          |      "attempt": 1,
                          |      "shardIndex": -1
                          |    }]
                          |  },
                          |  "id": "$workflowId"
                          |}""".stripMargin),
    ("query4", query4, s"""{
                          |  "key1": "value2",
                          |  "key2": "value1",
                          |  "calls": {
                          |    "sup.sup": [{
                          |      "dog": {},
                          |      "attempt": 1,
                          |      "shardIndex": -1
                          |    }]
                          |  },
                          |  "id": "$workflowId"
                          |}""".stripMargin),
    ("query5", query5, s"""{
                          |  "calls": {
                          |    "sup.sup": [{
                          |      "dog": {},
                          |      "attempt": 1,
                          |      "shardIndex": -1
                          |    }]
                          |  },
                          |  "id": "$workflowId"
                          |}""".stripMargin)
  )

  "MetadataServiceActor" should {

    testCases foreach { case (name, query, expectation) =>

      s"perform $name correctly" in {
        eventually(Timeout(10.seconds), Interval(2.seconds)) {
          val response = Await.result((actor ? GetMetadataAction(query)).mapTo[BuiltMetadataResponse], 1.seconds)

          response.response shouldBe expectation.parseJson
        }
      }

    }
  }

//      eventually(Timeout(10.seconds), Interval(2.seconds)) {
//        val response1 = Await.result((actor ? GetMetadataAction(query1)).mapTo[BuiltMetadataResponse], 1.seconds)
//        response1.response.prettyPrint shouldBe
//          s"""{
//             |  "${event1_2.key.key}": "${event1_2.value.map(_.value).get}",
//             |  "calls": {
//             |
//             |  },
//             |  "id": "$workflowId"
//             |}""".stripMargin
//      }
//
//      val response2 = Await.result((actor ? GetMetadataAction(query2)).mapTo[BuiltMetadataResponse], 10.seconds)
//      response2 shouldBe MetadataLookupResponse(query2, Seq(event2_1))
//
//      val response3 = Await.result((actor ? GetMetadataAction(query3)).mapTo[BuiltMetadataResponse], 10.seconds)
//      response3 shouldBe MetadataLookupResponse(query3, Seq(event3_1, event3_2))
//
//      val response4 = Await.result((actor ? GetMetadataAction(query4)).mapTo[BuiltMetadataResponse], 10.seconds)
//      response4 shouldBe MetadataLookupResponse(query4, Seq(event1_1, event1_2, event2_1, event3_1, event3_2))
//
//      val response5 = Await.result((actor ? GetMetadataAction(query5)).mapTo[BuiltMetadataResponse], 10.seconds)
//      response5 shouldBe MetadataLookupResponse(query5, Seq(event3_1, event3_2))

//      eventually(Timeout(10.seconds), Interval(2.seconds)) {
//      {
//        (for {
//          response1 <- (actor ? GetMetadataAction(query1)).mapTo[MetadataServiceResponse]
//          _ = response1 shouldBe MetadataLookupResponse(query1, Seq(event1_1, event1_2))
//
//          response2 <- (actor ? GetMetadataAction(query2)).mapTo[MetadataServiceResponse]
//          _ = response2 shouldBe MetadataLookupResponse(query2, Seq(event2_1))
//
//          response3 <- (actor ? GetMetadataAction(query3)).mapTo[MetadataServiceResponse]
//          _ = response3 shouldBe MetadataLookupResponse(query3, Seq(event3_1, event3_2))
//
//          response4 <- (actor ? GetMetadataAction(query4)).mapTo[MetadataServiceResponse]
//          _ = response4 shouldBe MetadataLookupResponse(query4, Seq(event1_1, event1_2, event2_1, event3_1, event3_2))
//
//          response5 <- (actor ? GetMetadataAction(query5)).mapTo[MetadataServiceResponse]
//          _ = response5 shouldBe MetadataLookupResponse(query5, Seq(event3_1, event3_2))
//
//        } yield ()).futureValue
//      }

}

object MetadataServiceActorSpec {
  val Config =
    """
      |services.MetadataService.db-batch-size = 3
      |services.MetadataService.db-flush-rate = 100 millis
    """.stripMargin
}
