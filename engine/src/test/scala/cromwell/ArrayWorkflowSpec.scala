package cromwell

import java.nio.file.{Files, Paths}
import java.util.UUID

import akka.testkit._
import wdl4s.WdlNamespaceWithWorkflow
import wdl4s.expression.NoFunctions
import wdl4s.types.{WdlArrayType, WdlFileType, WdlStringType}
import wdl4s.values.{WdlArray, WdlFile, WdlInteger, WdlString}
import cromwell.util.SampleWdl

import scala.language.postfixOps

class ArrayWorkflowSpec extends CromwellTestkitSpec {
  val tmpDir = Files.createTempDirectory("ArrayWorkflowSpec")
  val ns = WdlNamespaceWithWorkflow.load(SampleWdl.ArrayLiteral(tmpDir).wdlSource(""))
  val expectedArray = WdlArray(WdlArrayType(WdlFileType), Seq(WdlFile("f1"), WdlFile("f2"), WdlFile("f3")))

  "A task which contains a parameter " should {
    "accept an array for the value" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.ArrayIO,
        EventFilter.info(pattern = s"starting calls: wf.concat, wf.find, wf.serialize", occurrences = 1),
        expectedOutputs = Map(
          "wf.count_lines.count" -> WdlInteger(3),
          "wf.count_lines_array.count" -> WdlInteger(3),
          "wf.serialize.contents" -> WdlString("str1\nstr2\nstr3")
        )
      )
    }
  }

  "A static Array[File] declaration" should {
    "be a valid declaration" in {
      val declaration = ns.workflow.declarations.find(_.unqualifiedName == "arr").getOrElse {
        fail("Expected declaration 'arr' to be found")
      }
      val expression = declaration.expression.getOrElse {
        fail("Expected an expression for declaration 'arr'")
      }
      val value = expression.evaluate((s:String) => fail("No lookups"), NoFunctions).getOrElse {
        fail("Expected expression for 'arr' to evaluate")
      }
      value shouldEqual WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("f1"), WdlString("f2"), WdlString("f3")))
    }
    "be usable as an input" in {
      val catTask = ns.findTask("cat").getOrElse {
        fail("Expected to find task 'cat'")
      }
      val command = catTask.instantiateCommand(Map("files" -> expectedArray), NoFunctions).getOrElse {
        fail("Expected instantiation to work")
      }
      command shouldEqual "cat -s f1 f2 f3"
    }
    "Coerce Array[String] to Array[File] when running the workflow" in {
      val outputs = Map(
        "wf.cat.lines" -> WdlArray(WdlArrayType(WdlStringType), Seq(
            WdlString("line1"),
            WdlString("line2"),
            WdlString("line3"),
            WdlString("line4"),
            WdlString("line5")
          )
        )
      )
      val uuid = UUID.randomUUID()
      val sampleWdl = SampleWdl.ArrayLiteral(Paths.get("."))
      runWdlAndAssertOutputs(
        sampleWdl,
        eventFilter = EventFilter.info(pattern = s"starting calls: wf.cat", occurrences = 1),
        expectedOutputs = outputs
      )
      sampleWdl.cleanup()
    }
  }
}
