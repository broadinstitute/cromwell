package cromwell.engine.workflow.lifecycle.execution.callcaching

import cromwell.core.Tags.DbmsTest
import cromwell.core.WorkflowId
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql.joins.CallCachingJoin
import cromwell.database.sql.tables._
import cromwell.services.database.{DatabaseSystem, DatabaseTestKit, EngineDatabaseType}
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.time.{Millis, Seconds, Span}
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import org.specs2.mock.Mockito

import scala.concurrent.ExecutionContext

class CallCachingSlickDatabaseSpec
  extends FlatSpec with Matchers with ScalaFutures with BeforeAndAfterAll with Mockito with TableDrivenPropertyChecks {

  implicit val ec = ExecutionContext.global
  implicit val defaultPatience = PatienceConfig(scaled(Span(5, Seconds)), scaled(Span(100, Millis)))

  // Test with and without prefixes. With prefixing tests accessing the detritus value CLOB, especially with MariaDB.
  // https://jira.mariadb.org/browse/CONJ-717
  val allowResultReuseTests = Table(
    ("description", "prefixOption"),
    ("without prefixes", None),
    ("with some prefixes", Option(List("prefix1", "prefix2", "prefix3", "prefix4"))),
    ("with thousands of prefixes", Option((1 to 10000).map("prefix" + _).toList)),
  )

  DatabaseSystem.All foreach { databaseSystem =>
    behavior of s"CallCachingSlickDatabase on ${databaseSystem.shortName}"

    lazy val dataAccess = DatabaseTestKit.initializedDatabaseFromSystem(EngineDatabaseType, databaseSystem)

    forAll(allowResultReuseTests) { (description, prefixOption) =>

      val idA = WorkflowId.randomId().toString
      val callA = "AwesomeWorkflow.GoodJob"
      val callCachingEntryA = CallCachingEntry(
        idA,
        callA,
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

      val callCachingSimpletonsA = Seq(
        CallCachingSimpletonEntry("simpleKey", "simpleValue".toClobOption, "string")
      )

      val callCachingDetritusesA = Seq(
        CallCachingDetritusEntry("detritusKey", "detritusValue".toClobOption)
      )

      val aggregation = Option(CallCachingAggregationEntry("BASE_AGGREGATION", Option("FILE_AGGREGATION")))

      it should s"honor allowResultReuse $description" taggedAs DbmsTest in {
        (for {
          _ <- dataAccess.addCallCaching(Seq(
            CallCachingJoin(
              callCachingEntryA,
              callCachingHashEntriesA,
              aggregation,
              callCachingSimpletonsA, callCachingDetritusesA
            )
          ),
            100
          )
          hasBaseAggregation <- dataAccess.hasMatchingCallCachingEntriesForBaseAggregation(
            "BASE_AGGREGATION",
            prefixOption
          )
          _ = hasBaseAggregation shouldBe false
          hit <- dataAccess.findCacheHitForAggregation(
            "BASE_AGGREGATION",
            Option("FILE_AGGREGATION"),
            callCachePathPrefixes = prefixOption,
            1
          )
          _ = hit shouldBe empty
        } yield ()).futureValue
      }

      it should s"retrieve CallCacheJoin for call $description" taggedAs DbmsTest in {
        (for {
          join <- dataAccess.callCacheJoinForCall(idA, callA, 1)
          _ = join shouldBe defined
          getJoin = join.get
          // We can't compare directly because the ones out from the DB have IDs filled in, so just compare the relevant values
          _ = getJoin
            .callCachingHashEntries
            .map(e => (e.hashKey, e.hashValue)) should contain theSameElementsAs
            callCachingHashEntriesA.map(e => (e.hashKey, e.hashValue))
          _ = getJoin
            .callCachingSimpletonEntries
            .map(e => (e.simpletonKey, e.simpletonValue.map(_.toRawString))) should contain theSameElementsAs
            callCachingSimpletonsA.map(e => (e.simpletonKey, e.simpletonValue.map(_.toRawString)))
          _ = getJoin
            .callCachingAggregationEntry
            .map(e => (e.baseAggregation, e.inputFilesAggregation)) shouldBe
            aggregation.map(e => (e.baseAggregation, e.inputFilesAggregation))
          _ = getJoin
            .callCachingDetritusEntries
            .map(e => (e.detritusKey, e.detritusValue.map(_.toRawString))) should contain theSameElementsAs
            callCachingDetritusesA.map(e => (e.detritusKey, e.detritusValue.map(_.toRawString)))
        } yield ()).futureValue
      }

    }
  }
}
