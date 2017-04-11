package cromwell

import akka.testkit._
import cromwell.core.path.DefaultPathBuilder
import cromwell.util.SampleWdl
import wdl4s.expression.{NoFunctions, WdlFunctions}
import wdl4s.types.{WdlFileType, WdlIntegerType, WdlMapType, WdlStringType}
import wdl4s.values._
import wdl4s.{ImportResolver, WdlNamespaceWithWorkflow}

import scala.util.{Success, Try}

class MapWorkflowSpec extends CromwellTestKitWordSpec {
  private val pwd = DefaultPathBuilder.get(".")
  private val sampleWdl = SampleWdl.MapLiteral(pwd)
  val ns = WdlNamespaceWithWorkflow.load(sampleWdl.wdlSource(), Seq.empty[ImportResolver]).get
  val expectedMap = WdlMap(WdlMapType(WdlFileType, WdlStringType), Map(
    WdlFile("f1") -> WdlString("alice"),
    WdlFile("f2") -> WdlString("bob"),
    WdlFile("f3") -> WdlString("chuck")
  ))
  sampleWdl.cleanup()

  "A task which contains a parameter " should {
    "accept an array for the value" in {
      val sampleWdl = SampleWdl.MapLiteral(pwd)
      val callDir = "<<PWD>>/cromwell-executions/wf/<<UUID>>/call-write_map/inputs<<PWD>>"
      runWdlAndAssertOutputs(
        sampleWdl = sampleWdl,
        EventFilter.info(pattern = "Starting calls: wf.read_map:NA:1, wf.write_map:NA:1", occurrences = 1),
        expectedOutputs = Map(
          "wf.read_map.out_map" -> WdlMap(WdlMapType(WdlStringType, WdlIntegerType), Map(
            WdlString("x") -> WdlInteger(500),
            WdlString("y") -> WdlInteger(600),
            WdlString("z") -> WdlInteger(700)
          )),
          "wf.write_map.contents" -> WdlString(s"$callDir/f1\talice\n$callDir/f2\tbob\n$callDir/f3\tchuck")
        )
      )
      sampleWdl.cleanup()
    }
  }

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
      expectedMap.wdlType.coerceRawValue(value).get shouldEqual expectedMap
    }
    "be usable as an input" in {
      val writeMapTask = ns.findTask("write_map").getOrElse {
        fail("Expected to find task 'write_map'")
      }
      class CannedFunctions extends WdlFunctions[WdlValue] {
        def write_map(params: Seq[Try[WdlValue]]): Try[WdlFile] = Success(WdlFile("/test/map/path"))
        override def getFunction(name: String): WdlFunction = name match {
          case "write_map" => write_map
          case _ => throw new UnsupportedOperationException("Only write_map should be called")
        }
      }
      val command = writeMapTask.instantiateCommand(writeMapTask.inputsFromMap(Map("file_to_name" -> expectedMap)), new CannedFunctions).getOrElse {
        fail("Expected instantiation to work")
      }
      command shouldEqual "cat /test/map/path"
    }
    "Coerce Map[String, String] to Map[String, Int] when running the workflow" in {
      val sampleWdl = SampleWdl.MapLiteral(pwd)
      runWdlAndAssertOutputs(
        sampleWdl,
        eventFilter = EventFilter.info(pattern = "Starting calls: wf.read_map:NA:1, wf.write_map:NA:1", occurrences = 1),
        expectedOutputs = Map(
          "wf.read_map.out_map" -> WdlMap(WdlMapType(WdlStringType, WdlIntegerType), Map(
            WdlString("x") -> WdlInteger(500),
            WdlString("y") -> WdlInteger(600),
            WdlString("z") -> WdlInteger(700)
          ))
        )
      )
      sampleWdl.cleanup()
    }
  }
}
