package languages.wdl.draft2

import languages.wdl.draft2.MapWorkflowSpec._
import org.scalatest.{FlatSpec, Matchers}
import wdl.draft2.model.expression.{NoFunctions, WdlFunctions}
import wdl.draft2.model.{Draft2ImportResolver, WdlNamespaceWithWorkflow}
import wom.types.{WomMapType, WomSingleFileType, WomStringType}
import wom.values.{WomMap, WomSingleFile, WomString, WomValue}

import scala.util.{Success, Try}

class MapWorkflowSpec extends FlatSpec with Matchers {
  val namespace = WdlNamespaceWithWorkflow.load(WorkflowSource, Seq.empty[Draft2ImportResolver]).get
  val expectedMap = WomMap(WomMapType(WomSingleFileType, WomStringType), Map(
    WomSingleFile("f1") -> WomString("alice"),
    WomSingleFile("f2") -> WomString("bob"),
    WomSingleFile("f3") -> WomString("chuck")
  ))

  "A static Map[File, String] declaration" should "be a valid declaration" in {
    val declaration = namespace.workflow.declarations.find {
      _.unqualifiedName == "map"
    }.getOrElse {
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

  it should "be usable as an input" in {
    val writeMapTask = namespace.findTask("write_map").getOrElse {
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
}

object MapWorkflowSpec {
  val WorkflowSource =
    s"""
       |task write_map {
       |  Map[File, String] file_to_name
       |  command {
       |    cat $${write_map(file_to_name)}
       |  }
       |  output {
       |    String contents = read_string(stdout())
       |  }
       |}
       |
       |task read_map {
       |  command <<<
       |    python <<CODE
       |    map = {'x': 500, 'y': 600, 'z': 700}
       |    print("\\n".join(["{}\\t{}".format(k,v) for k,v in map.items()]))
       |    CODE
       |  >>>
       |  output {
       |    Map[String, Int] out_map = read_map(stdout())
       |  }
       |}
       |
       |workflow wf {
       |  Map[File, String] map = {"f1": "alice", "f2": "bob", "f3": "chuck"}
       |  call write_map {input: file_to_name = map}
       |  call read_map
       |}
      """.stripMargin
}
