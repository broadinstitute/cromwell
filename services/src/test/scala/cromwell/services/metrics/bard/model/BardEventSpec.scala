package cromwell.services.metrics.bard.model

import cromwell.services.metrics.bard.BardTestUtils
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class BardEventSpec extends AnyFlatSpec with Matchers with BardTestUtils {
  behavior of "BardEvent"

  it should "return a map of properties" in {
    val properties = taskSummaryEvent.getProperties
    properties.get("workflowId") should be(workflowId)
    properties.get("parentWorkflowId") should be(parentWorkflowId)
    properties.get("rootWorkflowId") should be(rootWorkflowId)
    properties.get("jobTag") should be(jobTag)
    properties.get("jobFullyQualifiedName") should be(jobFqn)
    properties.get("jobIndex") should be(jobIndex)
    properties.get("jobAttempt") should be(jobAttempt)
    properties.get("terminalState") should be(terminalState)
    properties.get("cloudPlatform") should be(cloudPlatform)
    properties.get("dockerImage") should be(dockerImage)
    properties.get("cpuCount") should be(cpu)
    properties.get("memoryBytes") should be(memory)
    properties.get("startTime") should be(start)
    properties.get("cpuStartTime") should be(cpuStart)
    properties.get("endTime") should be(end)
    properties.get("jobSeconds") should be(jobSeconds)
    properties.get("cpuSeconds") should be(cpuSeconds)
  }

  it should "return a map of properties with None properties set to null" in {
    val properties = taskSummaryEvent.copy(parentWorkflowId = None, jobIndex = None).getProperties
    properties.get("workflowId") should be(workflowId)
    Option(properties.get("parentWorkflowId")) should be(None)
    Option(properties.get("jobIndex")) should be(None)

    properties.containsKey("parentWorkflowId") should be(true)
    properties.containsKey("jobIndex") should be(true)
  }
}
