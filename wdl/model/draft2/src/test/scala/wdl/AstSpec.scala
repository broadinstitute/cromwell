package wdl

import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.draft2.model.{AstTools, WdlNamespace}

class AstSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {
  val namespace = WdlNamespace.loadUsingSource(SampleWdl.ThreeStep.workflowSource(), None, None).get

  "Parser" should "produce AST with 3 Task nodes" in {
    AstTools.findAsts(namespace.ast, "Task").size shouldEqual 3
  }

  it should "produce AST with 1 Workflow node" in {
    AstTools.findAsts(namespace.ast, "Workflow").size shouldEqual 1
  }

  it should "produce AST with 3 Call nodes in the Workflow node" in {
    AstTools
      .findAsts(
        AstTools.findAsts(namespace.ast, "Workflow").head,
        "Call"
      )
      .size shouldEqual 3
  }
}
