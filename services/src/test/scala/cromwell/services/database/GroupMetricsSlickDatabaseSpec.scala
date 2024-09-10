package cromwell.services.database

import com.dimafeng.testcontainers.Container
import cromwell.core.Tags.DbmsTest
import cromwell.database.sql.SqlConverters.OffsetDateTimeToSystemTimestamp
import cromwell.database.sql.tables.GroupMetricsEntry
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.time.{Millis, Seconds, Span}

import java.time.OffsetDateTime
import scala.concurrent.ExecutionContext

class GroupMetricsSlickDatabaseSpec extends AnyFlatSpec with Matchers with ScalaFutures {

  implicit val ec: ExecutionContext = ExecutionContext.global
  implicit val defaultPatience: PatienceConfig = PatienceConfig(scaled(Span(10, Seconds)), scaled(Span(100, Millis)))

  val testHogGroup1 = "groot-hog-group"
  val testHogGroup2 = "rocket-raccoon-hog-group"

  (DatabaseSystem.All diff List(HsqldbDatabaseSystem)) foreach { databaseSystem =>
    behavior of s"GroupMetricsSlickDatabase on ${databaseSystem.name}"

    val containerOpt: Option[Container] = DatabaseTestKit.getDatabaseTestContainer(databaseSystem)

    lazy val dataAccess =
      DatabaseTestKit.initializeDatabaseByContainerOptTypeAndSystem(containerOpt, EngineDatabaseType, databaseSystem)

    it should "start container if required" taggedAs DbmsTest in {
      containerOpt.foreach(_.start)
    }

    it should "add new row for group that doesn't exist in table" taggedAs DbmsTest in {
      (for {
        _ <- dataAccess.recordGroupMetricsEntry(GroupMetricsEntry(testHogGroup1, OffsetDateTime.now.toSystemTimestamp))
        rowCount <- dataAccess.countGroupMetricsEntries(testHogGroup1)
        _ = rowCount shouldBe 1
      } yield ()).futureValue
    }

    it should "update row for group that does exist in table" taggedAs DbmsTest in {
      // simulate group running into cloud quota delay for a long time
      (0 until 100).map { _ =>
        dataAccess.recordGroupMetricsEntry(GroupMetricsEntry(testHogGroup1, OffsetDateTime.now.toSystemTimestamp))
      }

      (for {
        rowCount <- dataAccess.countGroupMetricsEntries(testHogGroup1)
        // only 1 row for test hog group should exist in table
        _ = rowCount shouldBe 1
      } yield ()).futureValue
    }

    it should "record quota delay for 2 hog groups as expected in table" taggedAs DbmsTest in {
      (for {
        // groot-hog-group
        _ <- dataAccess.recordGroupMetricsEntry(GroupMetricsEntry(testHogGroup1, OffsetDateTime.now.toSystemTimestamp))
        // rocket-raccoon-hog-group
        _ <- dataAccess.recordGroupMetricsEntry(GroupMetricsEntry(testHogGroup2, OffsetDateTime.now.toSystemTimestamp))
        _ <- dataAccess.recordGroupMetricsEntry(GroupMetricsEntry(testHogGroup2, OffsetDateTime.now.toSystemTimestamp))
        // groot-hog-group
        _ <- dataAccess.recordGroupMetricsEntry(GroupMetricsEntry(testHogGroup1, OffsetDateTime.now.toSystemTimestamp))
        rowCountGroup1 <- dataAccess.countGroupMetricsEntries(testHogGroup1)
        rowCountGroup2 <- dataAccess.countGroupMetricsEntries(testHogGroup2)
        // only 1 row for each hog group should exist in table
        _ = rowCountGroup1 shouldBe 1
        _ = rowCountGroup2 shouldBe 1
      } yield ()).futureValue
    }

    it should "return correct group experiencing active quota exhaustion" taggedAs DbmsTest in {
      (for {
        // groot-hog-group records quota exhaustion at X-20 minutes
        _ <- dataAccess.recordGroupMetricsEntry(
          GroupMetricsEntry(testHogGroup1, OffsetDateTime.now.minusMinutes(20).toSystemTimestamp)
        )
        // rocket-raccoon-hog-group records quota exhaustion at X-18 minutes
        _ <- dataAccess.recordGroupMetricsEntry(
          GroupMetricsEntry(testHogGroup2, OffsetDateTime.now.minusMinutes(18).toSystemTimestamp)
        )
        // groot-hog-group records quota exhaustion again at X-10 minutes
        _ <- dataAccess.recordGroupMetricsEntry(GroupMetricsEntry(testHogGroup1, OffsetDateTime.now.toSystemTimestamp))
        // check that it returns only 'groot-hog-group' as currently experiencing quota exhaustion
        quotaExhaustedGroups <- dataAccess.getQuotaExhaustedGroups(
          OffsetDateTime.now.minusMinutes(15).toSystemTimestamp
        )
        _ = quotaExhaustedGroups.size shouldBe 1
        _ = quotaExhaustedGroups.toList.head shouldBe testHogGroup1
      } yield ()).futureValue
    }

    it should "return empty list if no group is experiencing active quota exhaustion" taggedAs DbmsTest in {
      (for {
        // groot-hog-group records quota exhaustion at X-20 minutes
        _ <- dataAccess.recordGroupMetricsEntry(
          GroupMetricsEntry(testHogGroup1, OffsetDateTime.now.minusMinutes(20).toSystemTimestamp)
        )
        // rocket-raccoon-hog-group records quota exhaustion at X-18 minutes
        _ <- dataAccess.recordGroupMetricsEntry(
          GroupMetricsEntry(testHogGroup2, OffsetDateTime.now.minusMinutes(18).toSystemTimestamp)
        )
        // check that it returns empty list
        quotaExhaustedGroups <- dataAccess.getQuotaExhaustedGroups(
          OffsetDateTime.now.minusMinutes(15).toSystemTimestamp
        )
        _ = quotaExhaustedGroups shouldBe empty
      } yield ()).futureValue
    }

    it should "close the database" taggedAs DbmsTest in {
      dataAccess.close()
    }

    it should "stop container if required" taggedAs DbmsTest in {
      containerOpt.foreach(_.stop())
    }
  }
}
