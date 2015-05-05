package cromwell

import cromwell.binding.WdlBinding
import org.scalatest.{Matchers, FlatSpec}

class AstSpec extends FlatSpec with Matchers with WdlThreeStepFixture {
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
