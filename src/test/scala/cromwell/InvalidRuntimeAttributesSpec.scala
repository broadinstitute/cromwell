package cromwell

import java.net.URL

import akka.testkit.{EventFilter, TestActorRef}
import cromwell.engine.ExecutionIndex._
import cromwell.engine._
import cromwell.engine.backend.jes.{JesAttributes, JesBackend}
import cromwell.engine.io.gcs.{GoogleConfiguration, Refresh, ServiceAccountMode, SimpleClientSecrets}
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

  "A workflow with a task with invalid runtime attributes" should {
    "fail on JES Backend" in {
      val jesBackend = new JesBackend(actorSystem) {
        private val anyString = ""
        private val anyURL: URL = null
        override lazy val jesConf = new JesAttributes(
          project = anyString,
          executionBucket = anyString,
          endpointUrl = anyURL) {
        }
        override def callGcsPath(descriptor: WorkflowDescriptor, callName: String, index: ExecutionIndex): String = "gs://fake/path"
        override def jesUserConnection(workflow: WorkflowDescriptor) = null
        override lazy val jesCromwellInterface = null
        override lazy val googleConf = GoogleConfiguration("appName", ServiceAccountMode("accountID", "pem"), None)
      }

      val workflowSources = WorkflowSourceFiles(SampleWdl.HelloWorld.wdlSource(), SampleWdl.HelloWorld.wdlJson, "{}")
      val submitMessage = WorkflowManagerActor.SubmitWorkflow(workflowSources)

      runWdlWithWorkflowManagerActor(
        wma = TestActorRef(new WorkflowManagerActor(jesBackend)),
        submitMsg = submitMessage,
        stdout = Map.empty,
        stderr = Map.empty,
        eventFilter = EventFilter.error(pattern = "RuntimeAttribute is not valid.", occurrences = 1),
        terminalState = WorkflowFailed
      )
    }

    "fail on Local Backend" in {
      runWdl(
        sampleWdl = SampleWdl.HelloWorld,
        runtime =
          """ runtime { wrongAttribute: "nop" }""".stripMargin,
        eventFilter = EventFilter.error(pattern = "RuntimeAttribute is not valid", occurrences = 1),
        terminalState = WorkflowFailed
      )
    }
  }

}
