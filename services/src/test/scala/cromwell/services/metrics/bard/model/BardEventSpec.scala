package cromwell.services.metrics.bard.model

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

import java.time.Instant
import java.time.temporal.ChronoUnit
import java.util.UUID

class BardEventSpec extends AnyFlatSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "BardEvent"

  it should "return a map of properties" in {
    val workflowId = UUID.randomUUID()
    val jobIdKey = "jobId"
    val taskState = "Complete"
    val bardEvent = TaskSummaryEvent(
      workflowId,
      jobIdKey,
      taskState,
      "gcp",
      "ubuntu",
      2,
      1024d,
      Instant.now().minus(1000, ChronoUnit.SECONDS).toString,
      Instant.now().minus(100, ChronoUnit.SECONDS).toString
    )

    val properties = bardEvent.getProperties
    properties.get("workflowId") shouldBe workflowId
    properties.get("jobIdKey") shouldBe jobIdKey
    properties.get("terminalState") shouldBe taskState
    properties.get("cloud") shouldBe "gcp"
    properties.get("dockerImage") shouldBe "ubuntu"
    properties.get("cpuCount") shouldBe 2
    properties.get("memoryBytes") shouldBe 1024d
    properties.get("startTime") shouldBe bardEvent.startTime
    properties.get("endTime") shouldBe bardEvent.endTime
    properties.get("pushToMixpanel") shouldBe false
    properties.get("event") shouldBe "task:summary"

  }
}
