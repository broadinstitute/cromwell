package cromwell

import akka.testkit._
import cromwell.core.Tags.PostWomTest
import cromwell.core.path.DefaultPathBuilder
import cromwell.util.SampleWdl
import wdl.expression.NoFunctions
import wdl.{ImportResolver, WdlNamespaceWithWorkflow}
import wom.types._
import wom.values._

class ArrayWorkflowSpec extends CromwellTestKitWordSpec {
  val tmpDir = DefaultPathBuilder.createTempDirectory("ArrayWorkflowSpec")
  val ns = WdlNamespaceWithWorkflow.load(SampleWdl.ArrayLiteral(tmpDir).workflowSource(), Seq.empty[ImportResolver]).get
  val expectedArray = WdlArray(WdlArrayType(WdlFileType), Seq(WdlFile("f1"), WdlFile("f2"), WdlFile("f3")))

  "A static Array[File] declaration" should {
    "be a valid declaration" in {
      val declaration = ns.workflow.declarations.find {_.unqualifiedName == "arr"}.getOrElse {
        fail("Expected declaration 'arr' to be found")
      }
      val expression = declaration.expression.getOrElse {
        fail("Expected an expression for declaration 'arr'")
      }
      val value = expression.evaluate((_: String) => fail("No lookups"), NoFunctions).getOrElse {
        fail("Expected expression for 'arr' to evaluate")
      }
      value shouldEqual WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("f1"), WdlString("f2"), WdlString("f3")))
    }
    "be usable as an input" in {
      val catTask = ns.findTask("cat").getOrElse {
        fail("Expected to find task 'cat'")
      }
      val command = catTask.instantiateCommand(catTask.inputsFromMap(Map("cat.files" -> expectedArray)), NoFunctions).getOrElse {
        fail("Expected instantiation to work")
      }
      command shouldEqual "cat -s f1 f2 f3"
    }
    // TODO WOM: should be fixed once call input evaluation is properly implemented
    "Coerce Array[String] to Array[File] when running the workflow" taggedAs PostWomTest ignore {
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
      val pwd = DefaultPathBuilder.get(".")
      val sampleWdl = SampleWdl.ArrayLiteral(pwd)
      runWdlAndAssertOutputs(
        sampleWdl,
        eventFilter = EventFilter.info(pattern = "Starting calls: wf.cat", occurrences = 1),
        expectedOutputs = outputs
      )
      sampleWdl.cleanup()
    }
  }
}
