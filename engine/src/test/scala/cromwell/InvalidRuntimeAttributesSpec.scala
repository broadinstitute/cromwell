package cromwell

import akka.testkit.EventFilter
import cromwell.core.WorkflowSucceeded
import cromwell.engine._
import cromwell.util.SampleWdl
import org.scalatest.BeforeAndAfterAll

class InvalidRuntimeAttributesSpec extends CromwellTestkitSpec with BeforeAndAfterAll {

  "A workflow with a task with one invalid runtime attribute" should {
    "succeed" in {
      runWdl(
        sampleWdl = SampleWdl.HelloWorld,
        runtime = """ runtime { wrongAttribute: "nop" }""".stripMargin,
        eventFilter = EventFilter.info(pattern = "transition from FinalizingWorkflowState to WorkflowSucceededState", occurrences = 1),
        terminalState = WorkflowSucceeded
      )
    }
  }

}