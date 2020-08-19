package cromwell.core.logging

import java.util.UUID

import cromwell.core.{PossiblyNotRootWorkflowId, RootWorkflowId}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class WorkflowLoggerSpec extends AnyFlatSpec with Matchers {
  private val workflowId = PossiblyNotRootWorkflowId(UUID.fromString("fc6cfad9-65e9-4eb7-853f-7e08c1c8cf8e"))
  private val rootWorkflowId = RootWorkflowId(UUID.fromString("04570ac7-39df-481c-8a37-eef6887f4a15"))

  behavior of "WorkflowLogger"

  it should "create a valid tag" in {
    val logger = new WorkflowLogger("LocalBackend", workflowId, rootWorkflowId, None)
    try {
      logger.tag should be("LocalBackend [UUID(fc6cfad9)]")
      logger.workflowLogPath shouldNot be(empty)
      logger.workflowLogPath.get.toString should
        endWith("cromwell-test-workflow-logs/workflow.04570ac7-39df-481c-8a37-eef6887f4a15.log")
    } finally {
      logger.close(true)
      ()
    }
  }
}
