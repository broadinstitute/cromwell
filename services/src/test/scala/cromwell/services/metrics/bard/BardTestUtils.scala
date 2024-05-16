package cromwell.services.metrics.bard

import cromwell.services.metrics.bard.model.TaskSummaryEvent

import java.time.Instant
import java.time.temporal.ChronoUnit
import java.util.UUID

trait BardTestUtils {
  val workflowId = UUID.randomUUID()
  val parentWorkflowId = UUID.randomUUID()
  val rootWorkflowId = UUID.randomUUID()
  val jobFqn = "workflowName.jobName"
  val jobIndex = 1
  val jobAttempt = 1
  val jobTag = s"$jobFqn:$jobIndex:$jobAttempt"
  val terminalState = "Complete"
  val cloud = "gcp"
  val dockerImage = "ubuntu"
  val cpu = 2
  val memory = 1024d
  val start = Instant.now().minus(1000, ChronoUnit.SECONDS).toString
  val end = Instant.now().minus(100, ChronoUnit.SECONDS).toString
  val taskSummaryEvent = TaskSummaryEvent(workflowId,
                                          Some(parentWorkflowId),
                                          rootWorkflowId,
                                          jobTag,
                                          jobFqn,
                                          Some(jobIndex),
                                          jobAttempt,
                                          terminalState,
                                          cloud,
                                          Some(dockerImage),
                                          cpu,
                                          memory,
                                          start,
                                          end
  )

}
