package wdl4s

import org.scalatest.{FlatSpec, Matchers}
import wdl4s.AstTools.AstNodeName
import wdl4s.parser.WdlParser.Ast

class ScopeSpec extends FlatSpec with Matchers {
  val namespace = WdlNamespaceWithWorkflow.load(SampleWdl.NestedScatterWdl.wdlSource()).get

  it should "throw an exception if trying to re-assign children on a scope" in {
    the [UnsupportedOperationException] thrownBy { namespace.workflow.children = Seq.empty } should have message "children is write-once"
  }

  it should "throw an exception if trying to generate a workflow from a non-workflow ast" in {
    val callAst: Ast = AstTools.findAsts(namespace.ast, AstNodeName.Call).head
    the [UnsupportedOperationException] thrownBy {
      Workflow(callAst, namespace.wdlSyntaxErrorFormatter)
    } should have message "Expecting Workflow AST, got a Call AST"
  }
}
