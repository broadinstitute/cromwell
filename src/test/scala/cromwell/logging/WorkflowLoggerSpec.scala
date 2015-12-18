package cromwell.logging

import java.util.UUID

import cromwell.CromwellTestkitSpec
import cromwell.binding.values.WdlValue
import cromwell.engine.backend.local.{LocalBackend, LocalBackendCall}
import cromwell.engine.workflow.CallKey
import cromwell.engine.{AbortRegistrationFunction, WorkflowDescriptor, WorkflowId, WorkflowSourceFiles}
import org.scalatest.{FlatSpec, Matchers}

class WorkflowLoggerSpec extends FlatSpec with Matchers {
  val testWorkflowManagerSystem = new CromwellTestkitSpec.TestWorkflowManagerSystem()

  val descriptor = WorkflowDescriptor(
    WorkflowId(UUID.fromString("fc6cfad9-65e9-4eb7-853f-7e08c1c8cf8e")),
    WorkflowSourceFiles(
      "task x {command {ps}} workflow w {call x}",
      "{}",
      "{}"
    )
  )
  val backend = LocalBackend(testWorkflowManagerSystem.actorSystem)
  val backendCall = LocalBackendCall(
    backend,
    descriptor,
    CallKey(descriptor.namespace.workflow.calls.find(_.unqualifiedName == "x").head, None),
    Map.empty[String, WdlValue],
    AbortRegistrationFunction(_ => ())
  )

  "WorkflowLogger" should "create a valid tag" in {
    backend.workflowLogger(descriptor).tag shouldBe "LocalBackend [UUID(fc6cfad9)]"
  }

  it should "create a valid tag for backend call" in {
    backend.workflowLoggerWithCall(backendCall).tag shouldBe "LocalBackend [UUID(fc6cfad9):x]"
  }
}
