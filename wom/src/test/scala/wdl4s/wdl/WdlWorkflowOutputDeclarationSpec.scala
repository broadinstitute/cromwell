package wdl4s.wdl

import org.scalatest.{Matchers, FlatSpec}

class WdlWorkflowOutputDeclarationSpec extends FlatSpec with Matchers {

  "WorkflowOutputDeclaration" should "match outputs" in {
    val declaration = WorkflowOutputWildcard("wf.mytask", wildcard = true, null)

    declaration.outputMatchesDeclaration("wf.mytask.a", wildcardsAllowed = true) shouldBe true
    declaration.outputMatchesDeclaration("wf.mytask.a", wildcardsAllowed = false) shouldBe false
    declaration.outputMatchesDeclaration("wf.mytaskwithadifferentname.a", wildcardsAllowed = true) shouldBe false
    declaration.outputMatchesDeclaration("wf.mytaskwithadifferentname.a", wildcardsAllowed = false) shouldBe false
  }

}
