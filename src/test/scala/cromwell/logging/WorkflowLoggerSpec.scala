package cromwell.logging

import java.util.UUID

import cromwell.CromwellTestkitSpec
import cromwell.engine.backend.BackendCallJobDescriptor
import wdl4s.values.WdlValue
import cromwell.engine.backend.local.{LocalBackend, LocalBackendCall}
import cromwell.engine.workflow.BackendCallKey
import cromwell.engine.{AbortRegistrationFunction, WorkflowDescriptor, WorkflowId, WorkflowSourceFiles}
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}

class WorkflowLoggerSpec extends FlatSpec with Matchers with BeforeAndAfterAll {
  val testWorkflowManagerSystem = new CromwellTestkitSpec.TestWorkflowManagerSystem()

  override protected def afterAll() = {
    testWorkflowManagerSystem.shutdownTestActorSystem()
    super.afterAll()
  }

  val descriptor = WorkflowDescriptor(
    WorkflowId(UUID.fromString("fc6cfad9-65e9-4eb7-853f-7e08c1c8cf8e")),
    WorkflowSourceFiles(
      "task x {command {ps}} workflow w {call x}",
      "{}",
      "{}"
    )
  )
  val backend = LocalBackend(testWorkflowManagerSystem.actorSystem)
  val jobDescriptor = BackendCallJobDescriptor(
    descriptor,
    BackendCallKey(descriptor.namespace.workflow.calls.find(_.unqualifiedName == "x").head, None, 1),
    Map.empty[String, WdlValue])

  val backendCall = LocalBackendCall(backend, jobDescriptor, callAbortRegistrationFunction = None)

  "WorkflowLogger" should "create a valid tag" in {
    backend.workflowLogger(descriptor).tag shouldBe "LocalBackend [UUID(fc6cfad9)]"
  }

  it should "create a valid tag for backend call" in {
    backend.workflowLoggerWithCall(backendCall).tag shouldBe "LocalBackend [UUID(fc6cfad9):x]"
  }
}
