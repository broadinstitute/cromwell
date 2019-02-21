package cromwell.services.metadata.impl

import java.sql.Timestamp
import java.util.UUID

import cromwell.database.sql.tables.{MetadataEntry, WorkflowMetadataSummaryEntry}
import org.scalatest.{FlatSpec, Matchers}
import javax.sql.rowset.serial.SerialClob
import common.assertion.ManyTimes._

import scala.util.Random

class MetadataSummarizationSpec extends FlatSpec with Matchers {

  behavior of "buildUpdatedSummary"

  def toClobOption(str: String) = Some(new SerialClob(str.toCharArray))

  it should "Correctly create an update" in {

    val workflowUuid = UUID.randomUUID().toString

    val existingSummary: Option[WorkflowMetadataSummaryEntry] = Some(WorkflowMetadataSummaryEntry(
      workflowExecutionUuid = workflowUuid,
      workflowName = None,
      workflowStatus = Some("Submitted"),
      startTimestamp = None,
      endTimestamp = None,
      submissionTimestamp = Some(Timestamp.valueOf("2019-02-21 18:59:51.777000"))
    ))

    val runningEvent = MetadataEntry(workflowUuid, None, None, None, "status", toClobOption("Running"), Some("string"), Timestamp.valueOf("2019-02-21 19:00:08.507000"))
    val workflowNameEvent = MetadataEntry(workflowUuid, None, None, None, "workflowName", toClobOption("test"), Some("string"), Timestamp.valueOf("2019-02-21 19:00:08.554000"))

    def makeEvents() = Seq(workflowNameEvent, runningEvent) ++ shuffleAndConcatenate(Seq(
      MetadataEntry(workflowUuid, None, None, None, "status", toClobOption("Submitted"), Some("string"), Timestamp.valueOf("2019-02-21 19:00:08.505000")),
      MetadataEntry(workflowUuid, None, None, None, "status", toClobOption("Starting"), Some("string"), Timestamp.valueOf("2019-02-21 19:00:08.505000")),
      runningEvent,
      workflowNameEvent,
    ), Seq.empty)

    val expectedSummary: WorkflowMetadataSummaryEntry = WorkflowMetadataSummaryEntry(
      workflowExecutionUuid = workflowUuid.toString,
      workflowName = Some("test"),
      workflowStatus = Some("Running"),
      startTimestamp = None,
      endTimestamp = None,
      submissionTimestamp = Some(Timestamp.valueOf("2019-02-21 18:59:51.777000"))
    )

    100.times {
      // Uses .toString() because object equality in timestamps seems a little off:
      val events = makeEvents()

      if (!(MetadataDatabaseAccess.buildUpdatedSummary(existingSummary, events).toString == expectedSummary.toString)) {
        println(s"Bad event sequence: ${events.map(MetadataSummarizationSpec.printable).mkString(System.lineSeparator, System.lineSeparator, System.lineSeparator)}")
        fail("Bad event sequence detected")
      }
    }
  }

  def shuffleAndConcatenate[A](as: Seq[A], result: Seq[A]): Seq[A] = {

    if (Random.nextInt(100) > 15) {
      val toInclude = as.filter(_ => Random.nextInt(100) > 40)
      val shuffled = Random.shuffle(result ++ toInclude)
      shuffleAndConcatenate(as, shuffled)
    } else {
      result
    }
  }

}

object MetadataSummarizationSpec {
  def printable(metadataEntry: MetadataEntry): String =
    s"""MetadataEntry("${metadataEntry.workflowExecutionUuid}", None, None, None, "${metadataEntry.metadataKey}", toClobOption("${metadataEntry.metadataValue.get.getSubString(1, metadataEntry.metadataValue.get.length.toInt)}"), Some("${metadataEntry.metadataValueType.get}"), Timestamp.valueOf("${metadataEntry.metadataTimestamp.toString}"))"""
}
