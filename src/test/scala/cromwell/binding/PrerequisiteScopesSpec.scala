package cromwell.binding

import cromwell.parser.BackendType
import cromwell.util.SampleWdl
import cromwell.util.SampleWdl.ScatterWdl
import org.scalatest.{FlatSpec, Matchers}

class PrerequisiteScopesSpec extends FlatSpec with Matchers {
  val namespace = NamespaceWithWorkflow.load((new ScatterWdl).wdlSource(), BackendType.LOCAL)
  val workflow = namespace.workflow
  val allCalls = workflow.collectAllCalls
  val allScatters = workflow.collectAllScatters

  def scopesByName(name: String) = workflow.callByName(name).get.prerequisiteScopes

  "ScatterWdl" should "have five calls" in {
    allCalls.size shouldEqual 5
  }

  it should "have one scatter block" in {
    allScatters.size shouldEqual 1
  }

  it should "have the scatter block with one prereq" in {
    allScatters.head.prerequisiteScopes.size shouldEqual 1
  }

  it should "have the scatter block with A as a prereq" in {
    val z = allScatters.head.prerequisiteScopes.head.fullyQualifiedName shouldEqual "w.A"
  }

  it should "have A not depend on anything" in {
    scopesByName("A") shouldBe empty
  }

  it should "have B depend on the scatter" in {
    val scopes = scopesByName("B")
    scopes.size shouldEqual 1
    scopes.head.name shouldBe "$scatter_0"
  }

  it should "have C depend on the scatter and B" in {
    val scopes = scopesByName("C")
    scopes.size shouldEqual 2
    scopes find { _.name == "$scatter_0" } shouldBe defined
    scopes find { _.name == "B"} shouldBe defined
  }

  it should "have D depend on B" in {
    val scopes = scopesByName("D")
    scopes.size shouldEqual 1
    scopes.head.name shouldBe "B"
  }

  it should "have E depend on the scatter" in {
    val scopes = scopesByName("E")
    scopes.size shouldEqual 1
    scopes.head.name shouldBe "$scatter_0"
  }
}
