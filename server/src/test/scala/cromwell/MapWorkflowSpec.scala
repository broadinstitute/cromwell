package cromwell

import akka.testkit._
import cromwell.core.path.DefaultPathBuilder
import cromwell.util.SampleWdl
import wdl.draft2.model.Draft2ImportResolver
import wdl.draft2.model.WdlNamespaceWithWorkflow
import wdl.draft2.model.expression.{NoFunctions, WdlFunctions}
import wom.types._
import wom.values._

import scala.util.{Success, Try}

class MapWorkflowSpec extends CromwellTestKitWordSpec {
  private val pwd = DefaultPathBuilder.get(".")
  private val sampleWdl = SampleWdl.MapLiteral(pwd)
  val ns = WdlNamespaceWithWorkflow.load(sampleWdl.workflowSource(), Seq.empty[Draft2ImportResolver]).get
  val expectedMap = WomMap(WomMapType(WomSingleFileType, WomStringType), Map(
    WomSingleFile("f1") -> WomString("alice"),
    WomSingleFile("f2") -> WomString("bob"),
    WomSingleFile("f3") -> WomString("chuck")
  ))
  sampleWdl.cleanup()

  "A static Map[File, String] declaration" should {
    "be a valid declaration" in {
      val declaration = ns.workflow.declarations.find {_.unqualifiedName == "map"}.getOrElse {
        fail("Expected declaration 'map' to be found")
      }
      val expression = declaration.expression.getOrElse {
        fail("Expected an expression for declaration 'map'")
      }
      val value = expression.evaluate((_: String) => fail("No lookups"), NoFunctions).getOrElse {
        fail("Expected expression for 'map' to evaluate")
      }
      expectedMap.womType.coerceRawValue(value).get shouldEqual expectedMap
    }
    "be usable as an input" in {
      val writeMapTask = ns.findTask("write_map").getOrElse {
        fail("Expected to find task 'write_map'")
      }
      class CannedFunctions extends WdlFunctions[WomValue] {
        def write_map(params: Seq[Try[WomValue]]): Try[WomSingleFile] = Success(WomSingleFile("/test/map/path"))
        override def getFunction(name: String): WdlFunction = name match {
          case "write_map" => write_map
          case _ => throw new UnsupportedOperationException("Only write_map should be called")
        }
      }
      val command = writeMapTask.instantiateCommand(writeMapTask.inputsFromMap(Map("file_to_name" -> expectedMap)), new CannedFunctions).getOrElse {
        fail("Expected instantiation to work")
      }
      command.head.commandString shouldEqual "cat /test/map/path"
    }
    "Coerce Map[String, String] to Map[String, Int] when running the workflow" in {
      val sampleWdl = SampleWdl.MapLiteral(pwd)
      runWdlAndAssertOutputs(
        sampleWdl,
        eventFilter = EventFilter.info(pattern = "Starting calls: wf.read_map:NA:1, wf.write_map:NA:1", occurrences = 1),
        expectedOutputs = Map(
          "wf.read_map.out_map" -> WomMap(WomMapType(WomStringType, WomIntegerType), Map(
            WomString("x") -> WomInteger(500),
            WomString("y") -> WomInteger(600),
            WomString("z") -> WomInteger(700)
          ))
        )
      )
      sampleWdl.cleanup()
    }
  }
}
