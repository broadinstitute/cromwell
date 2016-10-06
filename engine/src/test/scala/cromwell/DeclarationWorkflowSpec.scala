package cromwell

import wdl4s.types.{WdlFileType, WdlStringType}
import wdl4s.{NamespaceWithWorkflow, WorkflowInput}
import cromwell.util.SampleWdl
import org.scalatest.{Matchers, WordSpecLike}


class DeclarationWorkflowSpec extends Matchers with WordSpecLike {
  "A workflow with declarations in it" should {
    "compute inputs properly" in {
      NamespaceWithWorkflow.load(SampleWdl.DeclarationsWorkflow.wdlSource(runtime="")).workflow.inputs shouldEqual Map(
        "two_step.cat.file" -> WorkflowInput("two_step.cat.file", WdlFileType, postfixQuantifier = None),
        "two_step.cgrep.str_decl" -> WorkflowInput("two_step.cgrep.str_decl", WdlStringType, postfixQuantifier = None),
        "two_step.cgrep.pattern" -> WorkflowInput("two_step.cgrep.pattern", WdlStringType, postfixQuantifier = None),
        "two_step.flags_suffix" -> WorkflowInput("two_step.flags_suffix", WdlStringType, postfixQuantifier = None)
      )
    }
  }
}
