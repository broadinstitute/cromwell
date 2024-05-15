package cromwell.services.metrics.bard.model
import cromwell.core.WorkflowId
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

import java.time.Instant
import java.time.temporal.ChronoUnit
import java.util.UUID

class BardEventSpec extends AnyFlatSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "BardEvent"

  it should "return a map of properties" in {
    val workflowId = WorkflowId(UUID.randomUUID())
    val jobIdKey = "jobId"
    val taskState = "Complete"
    val bardEvent = TaskSummaryEvent(workflowId,
                                     jobIdKey,
                                     taskState,
                                     "gcp",
                                     "ubuntu",
                                     2,
                                     1024d,
                                     Instant.now().minus(1000, ChronoUnit.SECONDS),
                                     Instant.now().minus(100, ChronoUnit.SECONDS)
    )

    val properties = bardEvent.getProperties
    properties.get("pushToMixpanel") shouldBe false
    properties.get("event") shouldBe "task:summary"
    properties.get("terminationStatus") shouldBe "done"
    properties.get("returnCode") shouldBe 0

  }
}
