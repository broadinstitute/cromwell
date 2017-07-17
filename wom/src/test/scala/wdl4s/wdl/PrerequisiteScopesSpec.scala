package wdl4s.wdl

import wdl4s.wdl.SampleWdl.ScatterWdl
import org.scalatest.{FlatSpec, Matchers}

class PrerequisiteScopesSpec extends FlatSpec with Matchers {
  val namespace = WdlNamespaceWithWorkflow.load((new ScatterWdl).workflowSource(), Seq.empty).get
  val workflow = namespace.workflow
  val allCalls = workflow.calls
  val allScatters = workflow.scatters

  "ScatterWdl" should "have five calls" in {
    allCalls.size shouldEqual 5
  }

  it should "have one scatter block" in {
    allScatters.size shouldEqual 1
  }

  it should "have the scatter block with one prereq" in {
    allScatters.head.upstream.size shouldEqual 1
  }

  it should "have the scatter block with A as a prereq" in {
    allScatters.head.upstream.head.fullyQualifiedName shouldEqual "w.A"
  }

  it should "have A not depend on anything" in {
    workflow.calls.find(_.unqualifiedName == "A").get.upstream shouldBe empty
  }

  it should "have B.upstream == scatter" in {
    val callB = namespace.workflow.calls.find(_.unqualifiedName == "B").get
    callB.upstream shouldEqual Set(
      namespace.resolve("w.$scatter_0").get
    )
  }

  it should "have C.upstream == scatter, B" in {
    val callC = namespace.workflow.calls.find(_.unqualifiedName == "C").get
    callC.upstream shouldEqual Set(
      namespace.resolve("w.$scatter_0").get,
      namespace.resolve("w.B").get
    )
  }

  it should "have D.upstream == B" in {
    val callD = namespace.workflow.calls.find(_.unqualifiedName == "D").get
    callD.upstream shouldEqual Set(
      namespace.resolve("w.B").get
    )
  }

  it should "have E.upstream == scatter, A" in {
    val callE = namespace.workflow.calls.find(_.unqualifiedName == "E").get
    callE.upstream shouldEqual Set(
      namespace.resolve("w.$scatter_0").get
    )
  }
}
