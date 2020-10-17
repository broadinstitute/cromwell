package cromwell.core.logging

import java.util.UUID

import common.assertion.CromwellTimeoutSpec
import cromwell.core.{PossiblyNotRootWorkflowId, RootWorkflowId}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class JobLoggerSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {
  private val workflowId = PossiblyNotRootWorkflowId(UUID.fromString("fc6cfad9-65e9-4eb7-853f-7e08c1c8cf8e"))
  private val rootWorkflowId = RootWorkflowId(UUID.fromString("04570ac7-39df-481c-8a37-eef6887f4a15"))

  behavior of "JobLogger"

  it should "create a valid tag" in {
    val logger = new JobLogger("LocalBackend", workflowId, rootWorkflowId, "x", None)
    try {
      logger.tag should be("LocalBackend [UUID(fc6cfad9)x]")
      logger.workflowLogPath shouldNot be(empty)
      logger.workflowLogPath.get.toString should
        endWith("cromwell-test-workflow-logs/workflow.04570ac7-39df-481c-8a37-eef6887f4a15.log")
    } finally {
      logger.close(true)
      ()
    }
  }
}
