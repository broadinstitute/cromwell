package cromwell

import akka.testkit._
import cromwell.CromwellSpec.DockerTest
import cromwell.binding.types.{WdlStringType, WdlFileType}
import cromwell.binding.{WorkflowInput, NamespaceWithWorkflow}
import cromwell.binding.values.{WdlInteger, WdlString}
import cromwell.parser.BackendType
import cromwell.util.SampleWdl

import scala.language.postfixOps

class DeclarationWorkflowSpec extends CromwellTestkitSpec("DeclarationWorkflowSpec") {
  val outputs = Map(
    "two_step.cgrep.str" -> WdlString("foobar"),
    "two_step.cgrep.count" -> WdlInteger(1)
  )

  "A workflow with declarations in it" should {
    "compute inputs properly" in {
      NamespaceWithWorkflow.load(SampleWdl.DeclarationsWorkflow.wdlSource(runtime=""), BackendType.LOCAL).workflow.inputs shouldEqual Seq(
        WorkflowInput("two_step.cat.file", WdlFileType, postfixQuantifier = None),
        WorkflowInput("two_step.cgrep.pattern", WdlStringType, postfixQuantifier = None),
        WorkflowInput("two_step.cgrep.str_decl", WdlStringType, postfixQuantifier = None),
        WorkflowInput("two_step.flags_suffix", WdlStringType, postfixQuantifier = None)
      )
    }
    "honor the declarations" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.DeclarationsWorkflow,
        eventFilter = EventFilter.info(pattern = s"starting calls: two_step.cat", occurrences = 1),
        expectedOutputs = outputs
      )
    }
    "honor the declarations in a Docker environment" taggedAs DockerTest in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.DeclarationsWorkflow,
        eventFilter = EventFilter.info(pattern = s"starting calls: two_step.cat", occurrences = 1),
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
