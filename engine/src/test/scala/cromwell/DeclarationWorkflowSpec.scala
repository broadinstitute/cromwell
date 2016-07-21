package cromwell

import akka.testkit._
import cromwell.core.Tags.DockerTest
import wdl4s.types.{WdlStringType, WdlFileType}
import wdl4s.{WorkflowInput, NamespaceWithWorkflow}
import wdl4s.values.{WdlInteger, WdlString}
import cromwell.util.SampleWdl

import scala.language.postfixOps

class DeclarationWorkflowSpec extends CromwellTestkitSpec {
  val outputs = Map(
    "two_step.cgrep.str" -> WdlString("foobar"),
    "two_step.cgrep.count" -> WdlInteger(1)
  )

  "A workflow with declarations in it" should {
    "compute inputs properly" in {
      NamespaceWithWorkflow.load(SampleWdl.DeclarationsWorkflow.wdlSource(runtime="")).workflow.inputs shouldEqual Map(
        "two_step.cat.file" -> WorkflowInput("two_step.cat.file", WdlFileType, postfixQuantifier = None),
        "two_step.cgrep.str_decl" -> WorkflowInput("two_step.cgrep.str_decl", WdlStringType, postfixQuantifier = None),
        "two_step.cgrep.pattern" -> WorkflowInput("two_step.cgrep.pattern", WdlStringType, postfixQuantifier = None),
        "two_step.flags_suffix" -> WorkflowInput("two_step.flags_suffix", WdlStringType, postfixQuantifier = None)
      )
    }
    "honor the declarations" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.DeclarationsWorkflow,
        eventFilter = EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        expectedOutputs = outputs
      )
    }
    "honor the declarations in a Docker environment" taggedAs DockerTest in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.DeclarationsWorkflow,
        eventFilter = EventFilter.info(pattern = "Workflow complete", occurrences = 1),
        runtime =
          """runtime {
            |  docker: "ubuntu:latest"
            |}
          """.stripMargin,
        expectedOutputs = outputs
      )
    }
  }
}
