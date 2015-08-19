package cromwell.binding

import cromwell.binding.AstTools.AstNodeName
import cromwell.parser.BackendType
import cromwell.parser.WdlParser.Ast
import cromwell.util.SampleWdl
import org.scalatest.{FlatSpec, Matchers}

class ScopeSpec extends FlatSpec with Matchers {

  val namespace = NamespaceWithWorkflow.load(SampleWdl.ScatterWdl.wdlSource(), BackendType.LOCAL)

  it should "Have correct parent hierarchy" in {
    val calls: Seq[Call] = namespace.workflow.calls
    val scatters: Seq[Scatter] = namespace.workflow.scatters
    val scatter0: Scatter = scatters.find(_.name == "$scatter_0").get
    val scatter1: Scatter = scatters.find(_.name == "$scatter_1").get
    val scatter2: Scatter = scatters.find(_.name == "$scatter_2").get
    val scatter3: Scatter = scatters.find(_.name == "$scatter_3").get

    calls.find(_.fullyQualifiedName == "w.A").get.parent.get shouldEqual namespace.workflow
    calls.find(_.fullyQualifiedName == "w.B").get.parent.get shouldEqual scatter0
    calls.find(_.fullyQualifiedName == "w.C").get.parent.get shouldEqual scatter0
    calls.find(_.fullyQualifiedName == "w.E").get.parent.get shouldEqual scatter0
    calls.find(_.fullyQualifiedName == "w.G").get.parent.get shouldEqual scatter1
    calls.find(_.fullyQualifiedName == "w.H").get.parent.get shouldEqual scatter2
    calls.find(_.fullyQualifiedName == "w.F").get.parent.get shouldEqual scatter3
    calls.find(_.fullyQualifiedName == "w.D").get.parent.get shouldEqual namespace.workflow

    scatter0.parent.get shouldEqual namespace.workflow
    scatter1.parent.get shouldEqual scatter0
    scatter2.parent.get shouldEqual scatter0
    scatter3.parent.get shouldEqual namespace.workflow
    
    namespace.workflow.parent shouldEqual None
  }

  it should "Have correct children hierarchy" in {
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

  it should "throw an exception if trying to re-assign children on a scope" in {
    the [UnsupportedOperationException] thrownBy { namespace.workflow.setChildren(Seq.empty) } should have message "children is write-once"
  }

  it should "throw an exception if trying to generate a workflow from a non-workflow ast" in {
    val callAst: Ast = AstTools.findAsts(namespace.ast, AstNodeName.Call).head
    the [UnsupportedOperationException] thrownBy {
      Scope.generateWorkflow(callAst, namespace.namespaces, namespace.tasks, namespace.wdlSyntaxErrorFormatter)
    } should have message "Ast is not a 'Workflow Ast' but a 'Call Ast'"
  }

}
