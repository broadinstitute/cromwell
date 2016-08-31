package cromwell

import org.scalatest.{FlatSpec, Matchers, WordSpecLike}
import wdl4s.WdlNamespace
import wdl4s.expression.NoFunctions
import wdl4s.values.{WdlFile, WdlString}

import scala.language.postfixOps

class OptionalParamWorkflowSpec extends FlatSpec with Matchers with WordSpecLike {
  "A workflow with an optional parameter that has a prefix inside the tag" should {
    "not include that prefix if no value is specified" in {
      val wf = """
         |task find {
         |  String? pattern
         |  File root
         |  command {
         |    find ${root} ${"-name " + pattern}
         |  }
         |}
         |
         |workflow wf {
         |  call find
         |}
       """.stripMargin
      val ns = WdlNamespace.load(wf)
      val findTask = ns.findTask("find") getOrElse {
        fail("Expected to find task 'find'")
      }

      val instantiateWithoutValue = findTask.instantiateCommand(Map("root" -> WdlFile("src")), NoFunctions) getOrElse {
        fail("Expected instantiation to work")
      }
      instantiateWithoutValue shouldEqual "find src"

      val instantiateWithValue = findTask.instantiateCommand(Map(
        "root" -> WdlFile("src"),
        "pattern" -> WdlString("*.java")
      ), NoFunctions).getOrElse {fail("Expected instantiation to work")}
      instantiateWithValue shouldEqual "find src -name *.java"
    }
  }
}
