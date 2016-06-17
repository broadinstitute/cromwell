package cromwell.logging

import java.util.UUID

import cromwell.CromwellTestkitSpec
import cromwell.core.WorkflowId
import cromwell.engine.WorkflowSourceFiles
import cromwell.engine.backend.WorkflowDescriptorBuilder
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}
import wdl4s.values.WdlValue

class WorkflowLoggerSpec extends FlatSpec with Matchers with BeforeAndAfterAll with WorkflowDescriptorBuilder {
  lazy val testWorkflowManagerSystem = new CromwellTestkitSpec.TestWorkflowManagerSystem()
  override implicit lazy val actorSystem = testWorkflowManagerSystem.actorSystem

  override protected def afterAll() = {
    testWorkflowManagerSystem.shutdownTestActorSystem()
    super.afterAll()
  }

  // "WorkflowLogger" should "create a valid tag" in {
  ignore should "create a valid tag" in {
//    backend.workflowLogger(descriptor).tag shouldBe "LocalBackend [UUID(fc6cfad9)]"
  }

  ignore should "create a valid tag for backend call" in {
//    backend.jobLogger(jobDescriptor).tag shouldBe "LocalBackend [UUID(fc6cfad9):x]"
  }
}
