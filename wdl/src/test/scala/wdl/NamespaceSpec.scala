package wdl

import wom.values.WdlString

import scala.util.{Failure, Success, Try}

class NamespaceSpec extends WdlTest {

  "3 step WDL" should {
    val namespace = loadWdl("three_step/test.wdl")

    "coerceRawInputs (0)" in {
      namespace.coerceRawInputs(Map("three_step.cgrep.pattern" -> "abc")).get shouldEqual Map(
        "three_step.cgrep.pattern" -> WdlString("abc")
      )
    }

    "fail to coerceRawInputs if a required input is missing" in {
      namespace.coerceRawInputs(Map.empty[FullyQualifiedName, Any]) shouldBe a[Failure[_]]
    }
  }

  "WdlNamespace" should {
    "enforce optional output types" in {
      val namespace = Try(loadWdl("type_checks.wdl"))

      namespace match {
        case Failure(f) => f.getMessage should startWith("ERROR: oopsNotOptionalArray is declared as a Array[Int] but the expression evaluates to a Array[Int?]")
        case Success(_) => fail("Should have failed to load namespace")
      }
    }
  }
}
