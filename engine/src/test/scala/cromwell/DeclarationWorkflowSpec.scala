package cromwell

import wdl.{ImportResolver, WdlNamespaceWithWorkflow}
import cromwell.util.SampleWdl
import org.scalatest.{Matchers, WordSpecLike}
import wom.WorkflowInput
import wom.types.{WomFileType, WomStringType}


class DeclarationWorkflowSpec extends Matchers with WordSpecLike {
  "A workflow with declarations in it" should {
    "compute inputs properly" in {
      WdlNamespaceWithWorkflow.load(SampleWdl.DeclarationsWorkflow.workflowSource(runtime=""), Seq.empty[ImportResolver]).get.workflow.inputs shouldEqual Map(
        "two_step.cat.file" -> WorkflowInput("two_step.cat.file", WomFileType),
        "two_step.cgrep.str_decl" -> WorkflowInput("two_step.cgrep.str_decl", WomStringType),
        "two_step.cgrep.pattern" -> WorkflowInput("two_step.cgrep.pattern", WomStringType),
        "two_step.flags_suffix" -> WorkflowInput("two_step.flags_suffix", WomStringType)
      )
    }
  }
}
