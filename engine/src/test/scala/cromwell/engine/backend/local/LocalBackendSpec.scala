package cromwell.engine.backend.local

import akka.testkit.EventFilter
import cromwell.CromwellTestkitSpec
import cromwell.CromwellTestkitSpec.TestWorkflowManagerSystem
import cromwell.core.{WorkflowFailed, WorkflowSucceeded}
import cromwell.util.SampleWdl
import wdl4s.WdlSource

object LocalBackendSpec {
  object StdoutWdl extends SampleWdl {
    override def wdlSource(runtime: String): WdlSource =
      """task out { command { echo beach } RUNTIME }
        |workflow wf { call out }
      """.stripMargin.replaceAll("RUNTIME", runtime)
    override val rawInputs =  Map.empty[String, String]
  }

  object StderrWdl extends SampleWdl {
    override def wdlSource(runtime: String): WdlSource =
      """task err { command { echo beach >&2 } RUNTIME }
        |workflow wf { call err }
      """.stripMargin.replaceAll("RUNTIME", runtime)
    override val rawInputs =  Map.empty[String, String]
  }
}

class LocalBackendSpec extends CromwellTestkitSpec {
  import LocalBackendSpec._

  val testWorkflowManagerSystem = new TestWorkflowManagerSystem

  "LocalBackend" should {
    "allow stdout if failOnStderr is set" in {
      runWdl(
        sampleWdl = StdoutWdl,
        eventFilter = EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        runtime = "runtime {failOnStderr: true}",
        terminalState = WorkflowSucceeded)
    }

    "not allow stderr if failOnStderr is set" in {
      runWdl(
        sampleWdl = StderrWdl,
        eventFilter = EventFilter.info(pattern = "transitioning from FinalizingWorkflowState to WorkflowFailedState", occurrences = 1),
        runtime = "runtime {failOnStderr: true}",
        terminalState = WorkflowFailed)
    }

    "allow stdout if failOnStderr is not set" in {
      runWdl(
        sampleWdl = StdoutWdl,
        eventFilter = EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        runtime = "runtime {failOnStderr: false}",
        terminalState = WorkflowSucceeded)
    }

    "allow stderr if failOnStderr is not set" in {
      runWdl(
        sampleWdl = StderrWdl,
        eventFilter = EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        runtime = "runtime {failOnStderr: false}",
        terminalState = WorkflowSucceeded)
    }
  }
}
