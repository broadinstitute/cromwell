package cromwell

import org.scalatest.{Matchers, WordSpecLike}
import wdl.WdlNamespace
import wdl.expression.NoFunctions
import wom.values.{WdlFile, WdlString}


class OptionalParamWorkflowSpec extends Matchers with WordSpecLike {
  "A workflow with an optional parameter that has a prefix inside the tag" should {
    "not include that prefix if no value is specified" in {
      val wf = s"""
         |task find {
         |  String? pattern
         |  File root
         |  command {
         |    find $${root} $${"-name " + pattern}
         |  }
         |}
         |
         |workflow wf {
         |  call find
         |}
       """.stripMargin
      val ns = WdlNamespace.loadUsingSource(wf, None, None).get
      val findTask = ns.findTask("find") getOrElse {
        fail("Expected to find task 'find'")
      }

      val instantiateWithoutValue = findTask.instantiateCommand(findTask.inputsFromMap(Map("find.root" -> WdlFile("src"))), NoFunctions)
      instantiateWithoutValue.get shouldEqual "find src"

      val instantiateWithValue = findTask.instantiateCommand(findTask.inputsFromMap(Map(
        "find.root" -> WdlFile("src"),
        "find.pattern" -> WdlString("*.java")
      )), NoFunctions).getOrElse {fail("Expected instantiation to work")}
      instantiateWithValue shouldEqual "find src -name *.java"
    }
  }
}
