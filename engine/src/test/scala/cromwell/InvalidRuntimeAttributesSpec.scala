package cromwell

import akka.testkit.EventFilter
import cromwell.engine._
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
    "succeed" ignore {
      runWdl(
        sampleWdl = SampleWdl.HelloWorld,
        runtime = """ runtime { wrongAttribute: "nop" }""".stripMargin,
        eventFilter = EventFilter.info(pattern = "transitioning from Running to Succeeded", occurrences = 1),
        terminalState = WorkflowSucceeded
      )
    }
  }

}