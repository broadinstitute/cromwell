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

    it should "start container if required" taggedAs DbmsTest in {
      containerOpt.foreach(_.start)
    }

    forAll(allowResultReuseTests) { (description, prefixOption) =>
      val idA = WorkflowId.randomId().toString
      val callA = "AwesomeWorkflow.GoodJob"
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
                                             callCachingSimpletonsA,
                                             callCachingDetritusesA
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
            Set.empty,
            None
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
          _ = getJoin.callCachingHashEntries
            .map(e => (e.hashKey, e.hashValue)) should contain theSameElementsAs
            callCachingHashEntriesA.map(e => (e.hashKey, e.hashValue))
          _ = getJoin.callCachingSimpletonEntries
            .map(e => (e.simpletonKey, e.simpletonValue.map(_.toRawString))) should contain theSameElementsAs
            callCachingSimpletonsA.map(e => (e.simpletonKey, e.simpletonValue.map(_.toRawString)))
          _ = getJoin.callCachingAggregationEntry
            .map(e => (e.baseAggregation, e.inputFilesAggregation)) shouldBe
            aggregation.map(e => (e.baseAggregation, e.inputFilesAggregation))
          _ = getJoin.callCachingDetritusEntries
            .map(e => (e.detritusKey, e.detritusValue.map(_.toRawString))) should contain theSameElementsAs
            callCachingDetritusesA.map(e => (e.detritusKey, e.detritusValue.map(_.toRawString)))
        } yield ()).futureValue
      }
      it should s"not find cache hit when entry is older than maxResultAgeDays $description" taggedAs DbmsTest in {
        val idOld = WorkflowId.randomId().toString
        val callOld = "OldWorkflow.OldJob"

        // Create an entry that is 91 days old
        val oldTimestamp = Timestamp.from(Instant.now().minus(2, ChronoUnit.DAYS))
        val oldCallCachingEntry = CallCachingEntry(
          idOld,
          callOld,
          1,
          None,
          None,
          allowResultReuse = true,
          createdAt = oldTimestamp
        )

        val callCachingHashEntries = Seq(
          CallCachingHashEntry(
            hashKey = "input: String s1",
            hashValue = "HASH_OLD"
          )
        )

        val aggregationOld = Option(CallCachingAggregationEntry("BASE_AGG_OLD", Option("FILE_AGG_OLD")))

        (for {
          _ <- dataAccess.addCallCaching(
            Seq(
              CallCachingJoin(
                oldCallCachingEntry,
                callCachingHashEntries,
                aggregationOld,
                Seq.empty,
                Seq.empty
              )
            ),
            1
          )
          // Should find hit when maxResultAgeDays is None (no age filtering)
          hitWithoutAgeLimit <- dataAccess.findCacheHitForAggregation(
            "BASE_AGG_OLD",
            Option("FILE_AGG_OLD"),
            callCachePathPrefixes = prefixOption,
            Set.empty,
            maxResultAgeDays = None
          )
          _ = hitWithoutAgeLimit shouldBe defined

          // Should NOT find hit when maxResultAgeDays is 1 (entry is from 2 days ago)
          hitWithAgeLimit <- dataAccess.findCacheHitForAggregation(
            "BASE_AGG_OLD",
            Option("FILE_AGG_OLD"),
            callCachePathPrefixes = prefixOption,
            Set.empty,
            maxResultAgeDays = Some(1)
          )
          _ = hitWithAgeLimit shouldBe empty

          // Should find hit when maxResultAgeDays is 7 (entry is from 2 days ago)
          hitWithAgeLimit <- dataAccess.findCacheHitForAggregation(
            "BASE_AGG_OLD",
            Option("FILE_AGG_OLD"),
            callCachePathPrefixes = prefixOption,
            Set.empty,
            maxResultAgeDays = Some(7)
          )
          _ = hitWithAgeLimit shouldBe defined

        } yield ()).futureValue
      }

      it should s"return most recent cache entry first when multiple hits exist $description" taggedAs DbmsTest in {
        val baseAgg = "BASE_AGG_MULTIPLE"
        val fileAgg = Option("FILE_AGG_MULTIPLE")

        // Create three entries with different timestamps
        val now = Instant.now()
        val oldestEntry = CallCachingEntry(
          WorkflowId.randomId().toString,
          "Workflow.Job",
          1,
          None,
          None,
          allowResultReuse = true,
          createdAt = Timestamp.from(now.minus(3, ChronoUnit.DAYS))
        )

        val middleEntry = CallCachingEntry(
          WorkflowId.randomId().toString,
          "Workflow.Job",
          2,
          None,
          None,
          allowResultReuse = true,
          createdAt = Timestamp.from(now.minus(2, ChronoUnit.DAYS))
        )

        val newestEntry = CallCachingEntry(
          WorkflowId.randomId().toString,
          "Workflow.Job",
          3,
          None,
          None,
          allowResultReuse = true,
          createdAt = Timestamp.from(now.minus(1, ChronoUnit.DAYS))
        )

        val aggregation = Option(CallCachingAggregationEntry(baseAgg, fileAgg))

        (for {
          // Add entries in non-chronological order to verify sorting
          _ <- dataAccess.addCallCaching(
            Seq(
              CallCachingJoin(middleEntry, Seq.empty, aggregation, Seq.empty, Seq.empty),
              CallCachingJoin(oldestEntry, Seq.empty, aggregation, Seq.empty, Seq.empty),
              CallCachingJoin(newestEntry, Seq.empty, aggregation, Seq.empty, Seq.empty)
            ),
            3
          )

          // Should return the newest entry
          hit <- dataAccess.findCacheHitForAggregation(
            baseAgg,
            fileAgg,
            callCachePathPrefixes = prefixOption,
            Set.empty,
            maxResultAgeDays = None
          )
          _ = hit shouldBe defined

          // Verify it's the newest entry by checking the call cache entry ID
          newestJoin <- dataAccess.callCacheJoinForCall(
            newestEntry.workflowExecutionUuid,
            newestEntry.callFullyQualifiedName,
            newestEntry.jobIndex
          )
          _ = newestJoin shouldBe defined
          _ = hit.get shouldBe newestJoin.get.callCachingEntry.callCachingEntryId.get
        } yield ()).futureValue
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
