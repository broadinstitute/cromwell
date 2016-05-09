package cromwell.logging

import java.util.UUID

import cromwell.CromwellTestkitSpec
import cromwell.core.WorkflowId
import cromwell.engine.WorkflowSourceFiles
import cromwell.engine.backend.local.OldStyleLocalBackend
import cromwell.engine.backend.{OldStyleBackendCallJobDescriptor, WorkflowDescriptorBuilder}
import cromwell.engine.workflow.BackendCallKey
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import wdl4s.values.WdlValue

class WorkflowLoggerSpec extends FlatSpec with Matchers with BeforeAndAfterAll with WorkflowDescriptorBuilder {
  lazy val testWorkflowManagerSystem = new CromwellTestkitSpec.TestWorkflowManagerSystem()
  override implicit lazy val actorSystem = testWorkflowManagerSystem.actorSystem

  override protected def afterAll() = {
    testWorkflowManagerSystem.shutdownTestActorSystem()
    super.afterAll()
  }

  lazy val descriptor = materializeWorkflowDescriptorFromSources(id = WorkflowId(UUID.fromString("fc6cfad9-65e9-4eb7-853f-7e08c1c8cf8e")),
    workflowSources = WorkflowSourceFiles(
      "task x {command {ps}} workflow w {call x}",
      "{}",
      "{}"
    )
  )
  lazy val backend = OldStyleLocalBackend(CromwellTestkitSpec.DefaultLocalBackendConfigEntry, testWorkflowManagerSystem.actorSystem)
  lazy val jobDescriptor = OldStyleBackendCallJobDescriptor(
    descriptor,
    BackendCallKey(descriptor.namespace.workflow.calls.find(_.unqualifiedName == "x").head, None, 1),
    Map.empty[String, WdlValue])

  // "WorkflowLogger" should "create a valid tag" in {
  ignore should "create a valid tag" in {
    backend.workflowLogger(descriptor).tag shouldBe "LocalBackend [UUID(fc6cfad9)]"
  }

  ignore should "create a valid tag for backend call" in {
    backend.jobLogger(jobDescriptor).tag shouldBe "LocalBackend [UUID(fc6cfad9):x]"
  }
}
