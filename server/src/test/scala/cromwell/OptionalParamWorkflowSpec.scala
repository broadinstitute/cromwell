package cromwell

import common.validation.Validation._
import org.scalatest.matchers.should.Matchers
import org.scalatest.wordspec.AnyWordSpecLike
import wdl.draft2.model.WdlNamespace
import wdl.draft2.model.expression.NoFunctions
import wom.values.{WomSingleFile, WomString}

class OptionalParamWorkflowSpec extends Matchers with AnyWordSpecLike {
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

      val instantiateWithoutValue = findTask.instantiateCommand(
        findTask.inputsFromMap(Map("find.root" -> WomSingleFile("src"))),
        NoFunctions
      )
      instantiateWithoutValue.toTry.get.head.commandString shouldEqual "find src"

      val instantiateWithValue = findTask.instantiateCommand(findTask.inputsFromMap(Map(
        "find.root" -> WomSingleFile("src"),
        "find.pattern" -> WomString("*.java")
      )), NoFunctions).getOrElse {fail("Expected instantiation to work")}
      instantiateWithValue.head.commandString shouldEqual "find src -name *.java"
    }
  }
}
