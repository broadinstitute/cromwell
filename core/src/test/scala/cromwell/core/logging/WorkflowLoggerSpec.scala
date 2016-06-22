package cromwell.core.logging

import java.util.UUID

import cromwell.core.WorkflowId
import org.scalatest.{FlatSpec, Matchers}

class WorkflowLoggerSpec extends FlatSpec with Matchers {
  val workflowId = WorkflowId(UUID.fromString("fc6cfad9-65e9-4eb7-853f-7e08c1c8cf8e"))

   "WorkflowLogger" should "create a valid tag" in {
    new WorkflowLogger("LocalBackend", workflowId, None).tag shouldBe "LocalBackend [UUID(fc6cfad9)]"
  }
}
