package cromwell

import java.nio.file.{Files, Paths}

import akka.testkit._
import wdl4s.NamespaceWithWorkflow
import wdl4s.expression.{NoFunctions, WdlFunctions}
import wdl4s.types.{WdlFileType, WdlIntegerType, WdlMapType, WdlStringType}
import wdl4s.values._
import cromwell.util.SampleWdl

import scala.language.postfixOps
import scala.util.{Success, Try}

class MapWorkflowSpec extends CromwellTestkitSpec {
  val tmpDir = Files.createTempDirectory("MapWorkflowSpec")
  val ns = NamespaceWithWorkflow.load(SampleWdl.MapLiteral(Paths.get(".")).wdlSource(""))
  val expectedMap = WdlMap(WdlMapType(WdlFileType, WdlStringType), Map(
    WdlFile("f1") -> WdlString("alice"),
    WdlFile("f2") -> WdlString("bob"),
    WdlFile("f3") -> WdlString("chuck")
  ))

  "A task which contains a parameter " should {
    "accept an array for the value" in {
      runWdlAndAssertOutputs(
        sampleWdl = SampleWdl.MapLiteral(Paths.get(".")),
        EventFilter.info(pattern = s"starting calls: wf.read_map, wf.write_map", occurrences = 1),
        expectedOutputs = Map(
          "wf.read_map.out_map" -> WdlMap(WdlMapType(WdlStringType, WdlIntegerType), Map(
            WdlString("x") -> WdlInteger(500),
            WdlString("y") -> WdlInteger(600),
            WdlString("z") -> WdlInteger(700)
          )),
          "wf.write_map.contents" -> WdlString("f1\talice\nf2\tbob\nf3\tchuck")
        )
      )
    }
  }

  "A static Map[File, String] declaration" should {
    "be a valid declaration" in {
      val declaration = ns.workflow.declarations.find(_.unqualifiedName == "map").getOrElse {
        fail("Expected declaration 'map' to be found")
      }
      val expression = declaration.expression.getOrElse {
        fail("Expected an expression for declaration 'map'")
      }
      val value = expression.evaluate((s:String) => fail("No lookups"), NoFunctions).getOrElse {
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
      val command = writeMapTask.instantiateCommand(Map("file_to_name" -> expectedMap), new CannedFunctions).getOrElse {
        fail("Expected instantiation to work")
      }
      command shouldEqual "cat /test/map/path"
    }
    "Coerce Map[String, String] to Map[String, Int] when running the workflow" in {
      val outputs =
      runWdlAndAssertOutputs(
        SampleWdl.MapLiteral(Paths.get(".")),
        eventFilter = EventFilter.info(pattern = s"starting calls: wf.read_map, wf.write_map", occurrences = 1),
        expectedOutputs = Map(
          "wf.read_map.out_map" -> WdlMap(WdlMapType(WdlStringType, WdlIntegerType), Map(
            WdlString("x") -> WdlInteger(500),
            WdlString("y") -> WdlInteger(600),
            WdlString("z") -> WdlInteger(700)
          ))
        )
      )
    }
  }
}
