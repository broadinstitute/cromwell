package cromwell.binding

import org.scalatest.{FlatSpec, Matchers}

class AstSpec extends FlatSpec with Matchers with ThreeStepFixture {
  "Parser" should "produce AST with 3 Task nodes" in {
    WdlBinding.findAsts(binding.ast, "Task").size shouldEqual 3
  }

  it should "produce AST with 1 Workflow node" in {
    WdlBinding.findAsts(binding.ast, "Workflow").size shouldEqual 1
  }

  it should "produce AST with 3 Call nodes in the Workflow node" in {
    WdlBinding.findAsts(
      WdlBinding.findAsts(binding.ast, "Workflow").head,
      "Call"
    ).size shouldEqual 3
  }
}
