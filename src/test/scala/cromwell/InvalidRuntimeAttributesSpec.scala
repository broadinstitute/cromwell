package cromwell

import java.net.URL

import akka.testkit.{EventFilter, TestActorRef}
import cromwell.engine._
import cromwell.engine.backend.jes.{JesAttributes, JesBackend}
import cromwell.engine.io.gcs.{GoogleConfiguration, ServiceAccountMode}
import cromwell.engine.workflow.WorkflowManagerActor
import cromwell.util.SampleWdl
import org.scalatest.BeforeAndAfterAll

class InvalidRuntimeAttributesSpec extends CromwellTestkitSpec with BeforeAndAfterAll {

  val testWorkflowManagerSystem = new CromwellTestkitSpec.TestWorkflowManagerSystem()
  val actorSystem = testWorkflowManagerSystem.actorSystem

  override protected def afterAll() = {
    testWorkflowManagerSystem.shutdownTestActorSystem()
    super.afterAll()
  }

  "A workflow with a task with one invalid runtime attribute" should {
    "succeed" in {
      runWdl(
        sampleWdl = SampleWdl.HelloWorld,
        runtime =
          """ runtime { wrongAttribute: "nop" }""".stripMargin,
        eventFilter = EventFilter.info(pattern = "transitioning from Running to Succeeded", occurrences = 1),
        terminalState = WorkflowSucceeded
      )
    }
  }

  "A workflow with a task with invalid runtime attributes" should {
    "fail on JES Backend" in {
      val jesBackend = new JesBackend(actorSystem) {
        private val anyString = ""
        private val anyURL: URL = null
        override lazy val jesConf = new JesAttributes(
          project = anyString,
          executionBucket = anyString,
          endpointUrl = anyURL,
          maxPreemptionAttempts = 1,
          maxPollingInterval = 600) {
        }

        override def jesUserConnection(workflow: WorkflowDescriptor) = null

        override lazy val jesCromwellInterface = null
        override lazy val googleConf = GoogleConfiguration("appName", ServiceAccountMode("accountID", "pem"), None)
      }

      val workflowSources = WorkflowSourceFiles(SampleWdl.HelloWorld.wdlSource(), SampleWdl.HelloWorld.wdlJson, """ {"jes_gcs_root": "gs://fake/path"} """)

      runWdlWithWorkflowManagerActor(
        wma = TestActorRef(new WorkflowManagerActor(jesBackend)),
        sources = workflowSources,
        stdout = Map.empty,
        stderr = Map.empty,
        eventFilter = EventFilter.error(pattern = "RuntimeAttribute is not valid.", occurrences = 1),
        terminalState = WorkflowFailed,
        assertStdoutStderr = false
      )
    }
  }

}
