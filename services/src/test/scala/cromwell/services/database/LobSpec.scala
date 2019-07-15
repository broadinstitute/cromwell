package cromwell.services.database

import java.time.OffsetDateTime

import cromwell.core.Tags._
import cromwell.core.WorkflowId
import cromwell.database.sql.SqlConverters._
import cromwell.database.sql.joins.JobStoreJoin
import cromwell.database.sql.tables.{JobStoreEntry, JobStoreSimpletonEntry, WorkflowStoreEntry}
import javax.sql.rowset.serial.{SerialBlob, SerialException}
import org.scalatest.concurrent.ScalaFutures
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._

class LobSpec extends FlatSpec with Matchers with ScalaFutures {

  implicit val executionContext = ExecutionContext.global

  implicit val defaultPatience = PatienceConfig(timeout = scaled(5.seconds), interval = scaled(100.millis))

  DatabaseSystem.All foreach { databaseSystem =>

    behavior of s"CLOBs and BLOBs on ${databaseSystem.shortName}"

    lazy val database = DatabaseTestKit.initializedDatabaseFromSystem(EngineDatabaseType, databaseSystem)

    it should "fail to store and retrieve empty blobs" taggedAs DbmsTest in {
      // See notes in BytesToBlobOption
      import eu.timepit.refined.auto._
      val clob = "".toClob(default = "{}")
      val clobOption = "{}".toClobOption
      val emptyBlob = new SerialBlob(Array.empty[Byte])

      val workflowUuid = WorkflowId.randomId().toString
      val workflowStoreEntry = WorkflowStoreEntry(
        workflowExecutionUuid = workflowUuid,
        workflowRoot = Option("main"),
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        workflowDefinition = clobOption,
        workflowUrl = None,
        workflowInputs = clobOption,
        workflowOptions = clobOption,
        workflowState = "Submitted",
        cromwellId = None,
        heartbeatTimestamp = None,
        submissionTime = OffsetDateTime.now.toSystemTimestamp,
        importsZip = Option(emptyBlob),
        customLabels = clob,
        hogGroup = None,
      )

      val workflowStoreEntries = Seq(workflowStoreEntry)

      val future = databaseSystem match {
        case MysqlDatabaseSystem =>
          // MySQL crashes because it calls SerialBlob's getBytes instead of getBinaryStream
          database.addWorkflowStoreEntries(workflowStoreEntries).failed map { exception =>
            exception should be(a[SerialException])
            exception.getMessage should
              be("Invalid arguments: position cannot be less than 1 or greater than the length of the SerialBlob")
          }
        case _ =>
          database.addWorkflowStoreEntries(workflowStoreEntries)
      }

      future.futureValue
    }

    it should "store and retrieve empty clobs" taggedAs DbmsTest in {
      // See notes in StringToClobOption
      val workflowUuid = WorkflowId.randomId().toString
      val callFqn = "call.fqn"
      val jobIndex = 1
      val jobAttempt = 1
      val jobSuccessful = false
      val jobStoreEntry = JobStoreEntry(workflowUuid, callFqn, jobIndex, jobAttempt, jobSuccessful, None, None, None)
      val jobStoreSimpletonEntries = Seq(
        JobStoreSimpletonEntry("empty", "".toClobOption, "WdlString"),
        JobStoreSimpletonEntry("aEntry", "a".toClobOption, "WdlString")
      )
      val jobStoreJoins = Seq(JobStoreJoin(jobStoreEntry, jobStoreSimpletonEntries))

      val future = for {
        _ <- database.addJobStores(jobStoreJoins, 1)
        queried <- database.queryJobStores(workflowUuid, callFqn, jobIndex, jobAttempt)
        _ = {
          val jobStoreJoin = queried.get
          jobStoreJoin.jobStoreEntry.workflowExecutionUuid should be(workflowUuid)

          val emptyEntry = jobStoreJoin.jobStoreSimpletonEntries.find(_.simpletonKey == "empty").get
          emptyEntry.simpletonValue.toRawString should be("")

          val aEntry = jobStoreJoin.jobStoreSimpletonEntries.find(_.simpletonKey == "aEntry").get
          aEntry.simpletonValue.toRawString should be("a")
        }
        _ <- database.removeJobStores(Seq(workflowUuid))
      } yield ()
      future.futureValue
    }

    it should "store and retrieve empty blobs" taggedAs DbmsTest in {
      // See notes in BytesToBlobOption
      import eu.timepit.refined.auto._

      val submittedWorkflowState = "Submitted"
      val clob = "".toClob(default = "{}")
      val clobOption = "{}".toClobOption

      val emptyWorkflowUuid = WorkflowId.randomId().toString
      val emptyWorkflowStoreEntry = WorkflowStoreEntry(
        workflowExecutionUuid = emptyWorkflowUuid,
        workflowRoot = Option("main"),
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        workflowDefinition = clobOption,
        workflowUrl = None,
        workflowInputs = clobOption,
        workflowOptions = clobOption,
        workflowState = submittedWorkflowState,
        cromwellId = None,
        heartbeatTimestamp = None,
        submissionTime = OffsetDateTime.now.toSystemTimestamp,
        importsZip = Option(Array.empty[Byte]).toBlobOption,
        customLabels = clob,
        hogGroup = None,
      )

      val noneWorkflowUuid = WorkflowId.randomId().toString
      val noneWorkflowStoreEntry = WorkflowStoreEntry(
        workflowExecutionUuid = noneWorkflowUuid,
        workflowRoot = Option("main"),
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        workflowDefinition = clobOption,
        workflowUrl = None,
        workflowInputs = clobOption,
        workflowOptions = clobOption,
        workflowState = submittedWorkflowState,
        cromwellId = None,
        heartbeatTimestamp = None,
        submissionTime = OffsetDateTime.now.toSystemTimestamp,
        importsZip = None,
        customLabels = clob,
        hogGroup = None,
      )

      val aByte = 'a'.toByte
      val aByteWorkflowUuid = WorkflowId.randomId().toString
      val aByteWorkflowStoreEntry = WorkflowStoreEntry(
        workflowExecutionUuid = aByteWorkflowUuid,
        workflowRoot = Option("main"),
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        workflowDefinition = clobOption,
        workflowUrl = None,
        workflowInputs = clobOption,
        workflowOptions = clobOption,
        workflowState = submittedWorkflowState,
        cromwellId = None,
        heartbeatTimestamp = None,
        submissionTime = OffsetDateTime.now.toSystemTimestamp,
        importsZip = Option(Array(aByte)).toBlobOption,
        customLabels = clob,
        hogGroup = None,
      )

      val workflowStoreEntries = Seq(emptyWorkflowStoreEntry, noneWorkflowStoreEntry, aByteWorkflowStoreEntry)

      val future = for {
        _ <- database.addWorkflowStoreEntries(workflowStoreEntries)
        queried <- database.fetchWorkflowsInState(
          limit = Int.MaxValue,
          cromwellId = "crom-f00ba4",
          heartbeatTimestampTimedOut = 1.hour.ago,
          heartbeatTimestampTo = OffsetDateTime.now.toSystemTimestamp,
          workflowStateFrom = "Submitted",
          workflowStateTo = "Running",
          workflowStateExcluded = "On Hold"
        )

        _ = {
          val emptyEntry = queried.find(_.workflowExecutionUuid == emptyWorkflowUuid).get
          emptyEntry.importsZip.toBytesOption should be(None)

          val noneEntry = queried.find(_.workflowExecutionUuid == noneWorkflowUuid).get
          noneEntry.importsZip.toBytesOption should be(None)

          val aByteEntry = queried.find(_.workflowExecutionUuid == aByteWorkflowUuid).get
          aByteEntry.importsZip.toBytesOption.get.toSeq should be(Seq(aByte))
        }
        _ <- database.removeWorkflowStoreEntry(emptyWorkflowUuid)
        _ <- database.removeWorkflowStoreEntry(noneWorkflowUuid)
        _ <- database.removeWorkflowStoreEntry(aByteWorkflowUuid)
      } yield ()
      future.futureValue
    }

    it should "close the database" taggedAs DbmsTest in {
      database.close()
    }
  }
}
