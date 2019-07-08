package languages.wdl.draft2

import com.typesafe.config.ConfigFactory
import common.assertion.ErrorOrAssertions._
import common.validation.Validation._
import languages.wdl.draft2.ArrayCoercionsSpec._
import org.scalatest.{FlatSpec, Matchers}
import wdl.draft2.model.WdlNamespaceWithWorkflow
import wdl.draft2.model.expression.NoFunctions
import wom.core.WorkflowSource
import wom.expression.EmptyIoFunctionSet
import wom.types.{WomArrayType, WomSingleFileType, WomStringType}
import wom.values.{WomArray, WomSingleFile, WomString}


class ArrayCoercionsSpec extends FlatSpec with Matchers {

  var factory: WdlDraft2LanguageFactory = new WdlDraft2LanguageFactory(ConfigFactory.parseString(ConfigString))
  val arrayLiteralNamespace: WdlNamespaceWithWorkflow = WdlNamespaceWithWorkflow.load(ArrayDeclarationWorkflow, List.empty).get

  "A static Array[File] declaration" should "be a valid declaration" in {
    val declaration = arrayLiteralNamespace.workflow.declarations.find {_.unqualifiedName == "arr"}.getOrElse {
      fail("Expected declaration 'arr' to be found")
    }
    val expression = declaration.expression.getOrElse {
      fail("Expected an expression for declaration 'arr'")
    }
    expression.evaluate((_: String) => fail("No lookups"), NoFunctions).toChecked.shouldBeValid(
      WomArray(WomArrayType(WomStringType), Seq(WomString("f1"), WomString("f2"), WomString("f3")))
    )
  }

  "An Array[File]" should "be usable as an input" in {
    val expectedArray = WomArray(
      WomArrayType(WomSingleFileType),
      Seq(WomSingleFile("f1"), WomSingleFile("f2"), WomSingleFile("f3"))
    )
    val catTask = arrayLiteralNamespace.findTask("cat").getOrElse {
      fail("Expected to find task 'cat'")
    }
    val command = catTask.instantiateCommand(catTask.inputsFromMap(Map("cat.files" -> expectedArray)), NoFunctions).getOrElse {
      fail("Expected instantiation to work")
    }
    command.head.commandString shouldEqual "cat -s f1 f2 f3"
  }

  "An Array[String] input" should "coerce to an Array[File] input declaration" in {
    createExecutable(ArrayDeclarationWorkflow)
  }

  "An Array[Array[String]] input" should "coerce to Array[Array[File]] input declaration" in {
    createExecutable(ArrayOfArraysDeclarationWorkflow)
  }

  private def createExecutable(workflowSource: WorkflowSource): Unit = {
    val result = for {
      bundle <- factory.getWomBundle(workflowSource, None, "{}", List.empty, List.empty)
      _ <- factory.createExecutable(bundle, inputs = "{}", ioFunctions = new EmptyIoFunctionSet)
    } yield ()

    result match {
      case Left(errors) => fail("Bad things have happened: " + errors.toList.mkString(", "))
      case Right(_) =>
    }
  }
}

object ArrayCoercionsSpec {
  val ConfigString =
    """
      |{
      |  enabled: true
      |  strict-validation: true
      |}""".stripMargin

  val ArrayDeclarationWorkflow =
    s"""
       |task cat {
       |  Array[File]+ files
       |  command {
       |    cat -s $${sep = ' ' files}
       |  }
       |  output {
       |    Array[String] lines = read_lines(stdout())
       |  }
       |}
       |
       |workflow wf {
       |  Array[File] arr = ["f1", "f2", "f3"]
       |  call cat {input: files = arr}
       |}
      """.stripMargin

  val ArrayOfArraysDeclarationWorkflow =
    s"""
       |task subtask {
       |  Array[File] a
       |  command {
       |    cat $${sep = " " a}
       |  }
       |  output {
       |    String concatenated = read_string(stdout())
       |  }
       |}
       |
       |workflow wf {
       |  Array[Array[File]] nested_file = [["f1","f2"],["f3","f4"]]
       |
       |  scatter(n in nested_file) {
       |    call subtask {
       |      input: a = n
       |    }
       |  }
       |}
      """.stripMargin
}
