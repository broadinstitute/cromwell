package wdl

import scala.util.{Failure, Success, Try}

class NamespaceSpec extends WdlTest {

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
