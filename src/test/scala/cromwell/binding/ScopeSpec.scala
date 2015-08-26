package cromwell.binding

import cromwell.binding.AstTools.AstNodeName
import cromwell.parser.BackendType
import cromwell.parser.WdlParser.Ast
import cromwell.util.SampleWdl
import org.scalatest.{FlatSpec, Matchers}

class ScopeSpec extends FlatSpec with Matchers {

  val namespace = NamespaceWithWorkflow.load(SampleWdl.NestedScatterWdl.wdlSource(), BackendType.LOCAL)
  val calls: Seq[Call] = namespace.workflow.calls
  val scatters: Seq[Scatter] = namespace.workflow.scatters
  val scatter0: Scatter = scatters.find(_.name == "$scatter_0").get
  val scatter1: Scatter = scatters.find(_.name == "$scatter_1").get
  val scatter2: Scatter = scatters.find(_.name == "$scatter_2").get
  val scatter3: Scatter = scatters.find(_.name == "$scatter_3").get
  val callA: Call = calls.find(_.fullyQualifiedName == "w.A").get
  val callB: Call = calls.find(_.fullyQualifiedName == "w.B").get
  val callC: Call = calls.find(_.fullyQualifiedName == "w.C").get
  val callD: Call = calls.find(_.fullyQualifiedName == "w.D").get
  val callE: Call = calls.find(_.fullyQualifiedName == "w.E").get
  val callF: Call = calls.find(_.fullyQualifiedName == "w.F").get
  val callG: Call = calls.find(_.fullyQualifiedName == "w.G").get
  val callH: Call = calls.find(_.fullyQualifiedName == "w.H").get

  it should "Have correct parent hierarchy" in {
    callA.parent.get shouldEqual namespace.workflow
    callB.parent.get shouldEqual scatter0
    callC.parent.get shouldEqual scatter0
    callE.parent.get shouldEqual scatter0
    callG.parent.get shouldEqual scatter1
    callH.parent.get shouldEqual scatter2
    callF.parent.get shouldEqual scatter3
    callD.parent.get shouldEqual namespace.workflow

    scatter0.parent.get shouldEqual namespace.workflow
    scatter1.parent.get shouldEqual scatter0
    scatter2.parent.get shouldEqual scatter0
    scatter3.parent.get shouldEqual namespace.workflow

    namespace.workflow.parent shouldEqual None
  }

  it should "Have correct parent/child relationships" in {
    namespace.workflow.children shouldEqual Seq(callA, scatter0, scatter3, callD)

    scatter0.children shouldEqual Seq(callB, callC, callE, scatter1, scatter2)
    scatter1.children shouldEqual Seq(callG)
    scatter2.children shouldEqual Seq(callH)
    scatter3.children shouldEqual Seq(callF)

    callA.children shouldBe empty
    callB.children shouldBe empty
    callC.children shouldBe empty
    callD.children shouldBe empty
    callE.children shouldBe empty
    callF.children shouldBe empty
    callG.children shouldBe empty
    callH.children shouldBe empty
  }

  it should "Have correct ancestry for each Scope" in {
    callA.ancestry shouldEqual Seq(namespace.workflow)
    callB.ancestry shouldEqual Seq(scatter0, namespace.workflow)
    callC.ancestry shouldEqual Seq(scatter0, namespace.workflow)
    callE.ancestry shouldEqual Seq(scatter0, namespace.workflow)
    callG.ancestry shouldEqual Seq(scatter1, scatter0, namespace.workflow)
    callH.ancestry shouldEqual Seq(scatter2, scatter0, namespace.workflow)
    callF.ancestry shouldEqual Seq(scatter3, namespace.workflow)
    callD.ancestry shouldEqual Seq(namespace.workflow)
  }

  it should "Be able to determine common ancestor between two Scopes" in {
    callA.closestCommonAncestor(callH) shouldEqual Some(namespace.workflow)
    callH.closestCommonAncestor(callA) shouldEqual Some(namespace.workflow)
    callB.closestCommonAncestor(callC) shouldEqual Some(scatter0)
    callC.closestCommonAncestor(callB) shouldEqual Some(scatter0)
    callG.closestCommonAncestor(callH) shouldEqual Some(scatter0)
  }

  it should "throw an exception if trying to re-assign children on a scope" in {
    the [UnsupportedOperationException] thrownBy { namespace.workflow.children = Seq.empty } should have message "children is write-once"
  }

  it should "throw an exception if trying to generate a workflow from a non-workflow ast" in {
    val callAst: Ast = AstTools.findAsts(namespace.ast, AstNodeName.Call).head
    the [UnsupportedOperationException] thrownBy {
      Scope.generateWorkflow(callAst, namespace.namespaces, namespace.tasks, namespace.wdlSyntaxErrorFormatter)
    } should have message "Ast is not a 'Workflow Ast' but a 'Call Ast'"
  }
}
