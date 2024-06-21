package cromwell.services.metrics.bard

import akka.testkit.TestProbe
import cromwell.services.metrics.bard.model.TaskSummaryEvent

import java.time.OffsetDateTime
import java.time.temporal.ChronoUnit
import java.util.UUID

trait BardTestUtils {
  implicit val actorSystem = akka.actor.ActorSystem("BardTestUtils")
  val serviceRegistryProbe: TestProbe = TestProbe()
  val workflowId = UUID.randomUUID()
  val parentWorkflowId = UUID.randomUUID()
  val rootWorkflowId = UUID.randomUUID()
  val jobFqn = "workflowName.jobName"
  val jobIndex = 1
  val jobAttempt = 1
  val jobTag = s"$jobFqn:$jobIndex:$jobAttempt"
  val terminalState = "Complete"
  val platform = "gcp"
  val dockerImage = "ubuntu"
  val cpu = 2
  val memory = 1024d
  val start = OffsetDateTime.now().minus(1000, ChronoUnit.SECONDS).toString
  val cpuStart = OffsetDateTime.now().minus(800, ChronoUnit.SECONDS).toString
  val end = OffsetDateTime.now().minus(100, ChronoUnit.SECONDS).toString
  val jobSeconds = 900L
  val cpuSeconds = 700L
  val taskSummaryEvent = TaskSummaryEvent(
    workflowId,
    Some(parentWorkflowId),
    rootWorkflowId,
    jobTag,
    jobFqn,
    Some(jobIndex),
    jobAttempt,
    terminalState,
    Some(platform),
    Some(dockerImage),
    cpu,
    memory,
    start,
    Some(cpuStart),
    end,
    jobSeconds,
    Some(cpuSeconds)
  )

}
