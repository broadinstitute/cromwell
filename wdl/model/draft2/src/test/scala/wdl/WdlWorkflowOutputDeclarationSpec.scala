package wdl

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.draft2.model.WorkflowOutputWildcard

class WdlWorkflowOutputDeclarationSpec extends AnyFlatSpec with Matchers {

  "WorkflowOutputDeclaration" should "match outputs" in {
    val declaration = WorkflowOutputWildcard("wf.mytask", wildcard = true, null)

    declaration.outputMatchesDeclaration("wf.mytask.a", wildcardsAllowed = true) shouldBe true
    declaration.outputMatchesDeclaration("wf.mytask.a", wildcardsAllowed = false) shouldBe false
    declaration.outputMatchesDeclaration("wf.mytaskwithadifferentname.a", wildcardsAllowed = true) shouldBe false
    declaration.outputMatchesDeclaration("wf.mytaskwithadifferentname.a", wildcardsAllowed = false) shouldBe false
  }

}
