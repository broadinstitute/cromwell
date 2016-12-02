package cromwell

import wdl4s.types.{WdlFileType, WdlStringType}
import wdl4s.{ImportResolver, WdlNamespaceWithWorkflow, WorkflowInput}
import cromwell.util.SampleWdl
import org.scalatest.{Matchers, WordSpecLike}


class DeclarationWorkflowSpec extends Matchers with WordSpecLike {
  "A workflow with declarations in it" should {
    "compute inputs properly" in {
      WdlNamespaceWithWorkflow.load(SampleWdl.DeclarationsWorkflow.wdlSource(runtime=""), Seq.empty[ImportResolver]).workflow.inputs shouldEqual Map(
        "two_step.cat.file" -> WorkflowInput("two_step.cat.file", WdlFileType),
        "two_step.cgrep.str_decl" -> WorkflowInput("two_step.cgrep.str_decl", WdlStringType),
        "two_step.cgrep.pattern" -> WorkflowInput("two_step.cgrep.pattern", WdlStringType),
        "two_step.flags_suffix" -> WorkflowInput("two_step.flags_suffix", WdlStringType)
      )
    }
  }
}
