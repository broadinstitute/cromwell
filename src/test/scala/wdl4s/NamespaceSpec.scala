package wdl4s

import better.files._
import wdl4s.values.WdlString

import scala.util.Failure

class NamespaceSpec extends WdlTest {
  val threeStepWdl = "src/test/cases/three_step/test.wdl"

  threeStepWdl should {
    val namespace = loadWdlFile(File("src/test/cases/three_step/test.wdl"))

    "coerceRawInputs (0)" in {
      namespace.coerceRawInputs(Map("three_step.cgrep.pattern" -> "abc")).get shouldEqual Map(
        "three_step.cgrep.pattern" -> WdlString("abc")
      )
    }

    "fail to coerceRawInputs if a required input is missing" in {
      namespace.coerceRawInputs(Map.empty[FullyQualifiedName, Any]) shouldBe a[Failure[_]]
    }
  }
}
