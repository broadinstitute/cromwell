package cromwell.engine.workflow.lifecycle.execution.callcaching

import cats.data.NonEmptyList
import com.typesafe.config.ConfigFactory
import cromwell.core.Tags.DbmsTest
import cromwell.core.WorkflowId
import cromwell.database.slick.SlickDatabase
import cromwell.database.sql.joins.CallCachingJoin
import cromwell.database.sql.tables.{CallCachingAggregationEntry, CallCachingEntry, CallCachingHashEntry}
import cromwell.services.ServicesStore
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import org.specs2.mock.Mockito

import scala.concurrent.ExecutionContext

class CallCachingSlickDatabaseSpec extends FlatSpec with Matchers with ScalaFutures with BeforeAndAfterAll with Mockito {

  implicit val ec = ExecutionContext.global
  implicit val defaultPatience = PatienceConfig(scaled(Span(5, Seconds)), scaled(Span(100, Millis)))

  "SlickDatabase (hsqldb)" should behave like testWith("database")

  "SlickDatabase (mysql)" should behave like testWith("database-test-mysql")

  def testWith(configPath: String): Unit = {
    import ServicesStore.EnhancedSqlDatabase

    lazy val databaseConfig = ConfigFactory.load.getConfig(configPath)
    lazy val dataAccess = new SlickDatabase(databaseConfig).initialized

    val callCachingEntryA = CallCachingEntry(
      WorkflowId.randomId().toString,
      "AwesomeWorkflow.GoodJob",
      1,
      None,
      None,
      allowResultReuse = false
    )

    val callCachingHashEntriesA = Seq(
      CallCachingHashEntry(
        hashKey = "input: String s1",
        hashValue = "HASH_S1"
      ),
      CallCachingHashEntry(
        hashKey = "input: String s2",
        hashValue = "HASH_S2"
      ),
      CallCachingHashEntry(
        hashKey = "input: String s4",
        hashValue = "HASH_S4"
      )
    )

    val aggregation = Option(CallCachingAggregationEntry("BASE_AGGREGATION", Option("FILE_AGGREGATION")))

    it should "honor allowResultReuse" taggedAs DbmsTest in {
      (for {
        _ <- dataAccess.addCallCaching(Seq(
          CallCachingJoin(
            callCachingEntryA,
            callCachingHashEntriesA,
            aggregation,
            Seq.empty, Seq.empty
          )
        ),
          100
        )
        hasBaseAggregation <- dataAccess.hasMatchingCallCachingEntriesForBaseAggregation("BASE_AGGREGATION")
        _ = hasBaseAggregation shouldBe false
        hasHashPairMatch <- dataAccess.hasMatchingCallCachingEntriesForHashKeyValues(
          NonEmptyList.of("input: String s1" -> "HASH_S1")
        )
        _ = hasHashPairMatch shouldBe false
        hit <- dataAccess.findCacheHitForAggregation("BASE_AGGREGATION", Option("FILE_AGGREGATION"), 1)
        _ = hit shouldBe empty
      } yield ()).futureValue
    }
  }
}
