package wdl

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.draft2.model.AstTools.AstNodeName
import wdl.draft2.model.{AstTools, WdlNamespaceWithWorkflow, WdlWorkflow}
import wdl.draft2.parser.WdlParser.Ast

class ScopeSpec extends AnyFlatSpec with Matchers {
  val namespace = WdlNamespaceWithWorkflow.load(SampleWdl.NestedScatterWdl.workflowSource(), Seq.empty).get

  it should "throw an exception if trying to re-assign children on a scope" in {
    the [UnsupportedOperationException] thrownBy { namespace.workflow.children = Seq.empty } should have message "children is write-once"
  }

  it should "throw an exception if trying to generate a workflow from a non-workflow ast" in {
    val callAst: Ast = AstTools.findAsts(namespace.ast, AstNodeName.Call).head
    the [UnsupportedOperationException] thrownBy {
      WdlWorkflow(callAst, namespace.wdlSyntaxErrorFormatter)
    } should have message "Expecting Workflow AST, got a Call AST"
  }
}
