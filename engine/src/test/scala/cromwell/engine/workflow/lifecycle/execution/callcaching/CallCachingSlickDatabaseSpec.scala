package cromwell.engine.workflow.lifecycle.execution.callcaching

import com.dimafeng.testcontainers.Container
import common.assertion.CromwellTimeoutSpec
import cromwell.core.Tags.DbmsTest
import cromwell.core.WorkflowId
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql.joins.CallCachingJoin
import cromwell.database.sql.tables._
import cromwell.services.database.{DatabaseSystem, DatabaseTestKit, EngineDatabaseType}
import org.scalatest.BeforeAndAfterAll
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.time.{Millis, Seconds, Span}

import java.sql.Timestamp
import java.time.Instant
import java.time.temporal.ChronoUnit
import scala.concurrent.ExecutionContext

class CallCachingSlickDatabaseSpec
    extends AnyFlatSpec
    with CromwellTimeoutSpec
    with Matchers
    with ScalaFutures
    with BeforeAndAfterAll
    with TableDrivenPropertyChecks {

  implicit val ec: ExecutionContext = ExecutionContext.global
  implicit val defaultPatience: PatienceConfig = PatienceConfig(scaled(Span(5, Seconds)), scaled(Span(100, Millis)))

  // Test with and without prefixes. With prefixing tests accessing the detritus value CLOB, especially with MariaDB.
  // https://jira.mariadb.org/browse/CONJ-717
  private val allowResultReuseTests = Table(
    ("description", "prefixOption"),
    ("without prefixes", None),
    ("with some prefixes", Option(List("prefix1", "prefix2", "prefix3", "prefix4"))),
    ("with thousands of prefixes", Option((1 to 10000).map("prefix" + _).toList))
  )

  DatabaseSystem.All foreach { databaseSystem =>
    behavior of s"CallCachingSlickDatabase on ${databaseSystem.name}"

    val containerOpt: Option[Container] = DatabaseTestKit.getDatabaseTestContainer(databaseSystem)

    lazy val dataAccess =
      DatabaseTestKit.initializeDatabaseByContainerOptTypeAndSystem(containerOpt, EngineDatabaseType, databaseSystem)

    def createDetritusWithPrefix(prefix: Option[String], key: String, value: String): CallCachingDetritusEntry = {
      val fullValue = prefix match {
        case Some(p) => s"$p/$value"
        case None => value
      }
      val fullKey = prefix match {
        case Some(p) => s"$p/$key"
        case None => key
      }
      CallCachingDetritusEntry(
        detritusKey = fullKey,
        detritusValue = fullValue.toClobOption
      )
    }

    def createDetritusesBasedOnPrefixes(prefixOption: Option[List[String]]): Seq[CallCachingDetritusEntry] =
      prefixOption match {
        case Some(prefixes) =>
          prefixes.map { prefix =>
            createDetritusWithPrefix(Some(prefix), "simpleKey", "simpleValue")
          }
        case None =>
          Seq(createDetritusWithPrefix(None, "simpleKey", "simpleValue"))
      }

    it should "start container if required" taggedAs DbmsTest in {
      containerOpt.foreach(_.start)
    }

    forAll(allowResultReuseTests) { (description, prefixOption) =>
      // Create test data - the same test data is utilized for all tests below, but is unique to each (databaseSystem, prefixOption)

      // Entry A - allowResultReuse = false
      val idA = WorkflowId.randomId().toString
      val callA = "AwesomeWorkflow.GoodJobA"
      val callCachingEntryA = CallCachingEntry(
        idA,
        callA,
        1,
        None,
        None,
        allowResultReuse = false,
        createdAt = Timestamp.from(Instant.now())
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
      val aggregationA = Option(CallCachingAggregationEntry("BASE_AGGREGATION_A", Option("FILE_AGGREGATION_A")))

      // Entry B - allowResultReuse = true
      val idB = WorkflowId.randomId().toString
      val callB = "AwesomeWorkflow.GoodJobB"
      val callCachingEntryB = CallCachingEntry(
        idB,
        callB,
        1,
        None,
        None,
        allowResultReuse = true,
        createdAt = Timestamp.from(Instant.now())
      )

      val callCachingHashEntriesB = Seq(
        CallCachingHashEntry(
          hashKey = "input: String s5",
          hashValue = "HASH_S5"
        )
      )
      val aggregationB = Option(CallCachingAggregationEntry("BASE_AGGREGATION_B", Option("FILE_AGGREGATION_B")))

      // same across A and B
      val callCachingSimpletons = Seq(
        CallCachingSimpletonEntry("simpleKey", "simpleValue".toClobOption, "string")
      )
      val callCachingDetrituses = createDetritusesBasedOnPrefixes(prefixOption)

      // Create older entries
      val callOld = "OldWorkflow.OldJob"

      // Entry C - 1 day old
      val callCachingEntryC = CallCachingEntry(
        WorkflowId.randomId().toString,
        callOld,
        1,
        None,
        None,
        allowResultReuse = true,
        createdAt = Timestamp.from(Instant.now().minus(1, ChronoUnit.DAYS))
      )
      val callCachingHashEntriesC = Seq(
        CallCachingHashEntry(
          hashKey = "input: String s6",
          hashValue = "HASH_S6"
        )
      )

      // Entry D - 2 days old
      val callCachingEntryD = CallCachingEntry(
        WorkflowId.randomId().toString,
        callOld,
        2,
        None,
        None,
        allowResultReuse = true,
        createdAt = Timestamp.from(Instant.now().minus(2, ChronoUnit.DAYS))
      )
      val callCachingHashEntriesD = Seq(
        CallCachingHashEntry(
          hashKey = "input: String s7",
          hashValue = "HASH_S7"
        )
      )

      // same across old entries
      val aggregationOld = Option(CallCachingAggregationEntry("AGG_OLD", Option("FILE_AGG_OLD")))
      val callCachingDetritusesOld = createDetritusesBasedOnPrefixes(prefixOption)

      println(s"Prefix: $prefixOption")
      println(s"CallCachingDetrituses: $callCachingDetrituses")
      println(s"CallCachingDetritusesOld: $callCachingDetritusesOld")

      it should s"seed the database with test data $description" taggedAs DbmsTest in {
        // add call caching entries to DB with allowResultReuse = false
        dataAccess
          .addCallCaching(
            Seq(
              // with allowResultReuse = false
              CallCachingJoin(callCachingEntryA,
                              callCachingHashEntriesA,
                              aggregationA,
                              callCachingSimpletons,
                              callCachingDetrituses
              ),
              // with allowResultReuse = true
              CallCachingJoin(callCachingEntryB,
                              callCachingHashEntriesB,
                              aggregationB,
                              callCachingSimpletons,
                              callCachingDetrituses
              ),
              // one day old
              CallCachingJoin(callCachingEntryC,
                              callCachingHashEntriesC,
                              aggregationOld,
                              Seq.empty,
                              callCachingDetritusesOld
              ),
              // two days old
              CallCachingJoin(callCachingEntryD,
                              callCachingHashEntriesD,
                              aggregationOld,
                              Seq.empty,
                              callCachingDetritusesOld
              )
            ),
            100
          )
          .futureValue
      }

      it should s"honor allowResultReuse $description" taggedAs DbmsTest in {
        (for {
          hasBaseAggregation <- dataAccess.hasMatchingCallCachingEntriesForBaseAggregation(
            "BASE_AGGREGATION_A",
            prefixOption
          )
          _ = hasBaseAggregation shouldBe false
          hit <- dataAccess.findCacheHitForAggregation(
            "BASE_AGGREGATION_A",
            Option("FILE_AGGREGATION_A"),
            callCachePathPrefixes = prefixOption,
            Set.empty,
            None
          )
          _ = hit shouldBe empty
        } yield ()).futureValue
      }

      it should s"find a callCaching entry if allowResultReuse is true $description" taggedAs DbmsTest in {
        (for {
          hasBaseAggregation <- dataAccess.hasMatchingCallCachingEntriesForBaseAggregation(
            "BASE_AGGREGATION_B",
            prefixOption
          )
          _ = hasBaseAggregation shouldBe true
          hit <- dataAccess.findCacheHitForAggregation(
            "BASE_AGGREGATION_B",
            Option("FILE_AGGREGATION_B"),
            None,
            Set.empty,
            None
          )
          _ = hit shouldBe defined
        } yield ()).futureValue
      }

      it should s"retrieve CallCacheJoin for call $description" taggedAs DbmsTest in {
        (for {
          join <- dataAccess.callCacheJoinForCall(idA, callA, 1)
          _ = join shouldBe defined
          getJoin = join.get
          // We can't compare directly because the ones out from the DB have IDs filled in, so just compare the relevant values
          _ = getJoin.callCachingHashEntries
            .map(e => (e.hashKey, e.hashValue)) should contain theSameElementsAs
            callCachingHashEntriesA.map(e => (e.hashKey, e.hashValue))
          _ = getJoin.callCachingSimpletonEntries
            .map(e => (e.simpletonKey, e.simpletonValue.map(_.toRawString))) should contain theSameElementsAs
            callCachingSimpletons.map(e => (e.simpletonKey, e.simpletonValue.map(_.toRawString)))
          _ = getJoin.callCachingAggregationEntry
            .map(e => (e.baseAggregation, e.inputFilesAggregation)) shouldBe
            aggregationA.map(e => (e.baseAggregation, e.inputFilesAggregation))
          _ = getJoin.callCachingDetritusEntries
            .map(e => (e.detritusKey, e.detritusValue.map(_.toRawString))) should contain theSameElementsAs
            callCachingDetrituses.map(e => (e.detritusKey, e.detritusValue.map(_.toRawString)))
        } yield ()).futureValue
      }

      it should s"filter cache cache entry results based on maxResultAgeDays $description" taggedAs DbmsTest in {
        (for {
          // Should find hit when maxResultAgeDays is None (no age filtering)
          hitWithoutAgeLimit <- dataAccess.findCacheHitForAggregation(
            "AGG_OLD",
            Option("FILE_AGG_OLD"),
            None,
            Set.empty,
            maxResultAgeDays = None
          )
          _ = hitWithoutAgeLimit shouldBe defined

          // Should NOT find hit when maxResultAgeDays is 1 (entry is from 2 days ago)
          hitWithAgeLimit <- dataAccess.findCacheHitForAggregation(
            "AGG_OLD",
            Option("FILE_AGG_OLD"),
            None,
            Set.empty,
            maxResultAgeDays = Some(1)
          )
          _ = hitWithAgeLimit shouldBe empty

          // Should find hit when maxResultAgeDays is 7 (entry is from 2 days ago)
          hitWithAgeLimit <- dataAccess.findCacheHitForAggregation(
            "AGG_OLD",
            Option("FILE_AGG_OLD"),
            None,
            Set.empty,
            maxResultAgeDays = Some(7)
          )
          _ = hitWithAgeLimit shouldBe defined

        } yield ()).futureValue
      }

      it should s"return most recent cache entry first when multiple hits exist $description" taggedAs DbmsTest in {

        (for {
          // Should return the newest entry
          hit <- dataAccess.findCacheHitForAggregation(
            "AGG_OLD",
            Option("FILE_AGG_OLD"),
            None,
            Set.empty,
            maxResultAgeDays = None
          )

          _ = info(s"hit: ${hit.toString}")

          _ = hit shouldBe defined

          // Verify it's the newest entry by checking the call cache entry ID
          newestJoin <- dataAccess.callCacheJoinForCall(
            callCachingEntryD.workflowExecutionUuid,
            callCachingEntryD.callFullyQualifiedName,
            callCachingEntryD.jobIndex
          )
          _ = newestJoin shouldBe defined
          _ = hit.get shouldBe newestJoin.get.callCachingEntry.callCachingEntryId.get
        } yield ()).futureValue
      }

      it should s"clear entries from the database $description" taggedAs DbmsTest in {
        import dataAccess.dataAccess.driver.api._

        dataAccess.database
          .run(
            DBIO
              .seq(
                dataAccess.dataAccess.callCachingDetritusEntries.delete,
                dataAccess.dataAccess.callCachingSimpletonEntries.delete,
                dataAccess.dataAccess.callCachingHashEntries.delete,
                dataAccess.dataAccess.callCachingAggregationEntries.delete
              )
              .transactionally
          )
          .futureValue
      }
    }

    it should "close the database" taggedAs DbmsTest in {
      dataAccess.close()
    }

    it should "stop container if required" taggedAs DbmsTest in {
      containerOpt.foreach(_.stop())
    }
  }
}
