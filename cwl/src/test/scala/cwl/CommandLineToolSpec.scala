package cwl

import cats.instances.list._
import cats.syntax.traverse._
import common.validation.ErrorOr.ErrorOr
import org.scalatest.{FlatSpec, Matchers, ParallelTestExecution}
import wom.callable.Callable.RequiredInputDefinition
import wom.callable.RuntimeEnvironment
import wom.expression.IoFunctionSet
import wom.types._
import wom.values.{WomArray, WomBoolean, WomEvaluatedCallInputs, WomFile, WomFloat, WomInteger, WomObject, WomSingleFile, WomString, WomValue}

import scala.concurrent.Future
import scala.util.Try

class CommandLineToolSpec extends FlatSpec with Matchers with ParallelTestExecution {

  behavior of "CommandLineTool"

  val inputs: WomEvaluatedCallInputs = Map(
    RequiredInputDefinition("a", WomStringType) -> WomString("helloA"),
    RequiredInputDefinition("b", WomStringType) -> WomString("helloB"),
    RequiredInputDefinition("c", WomStringType) -> WomString("helloC"),
    RequiredInputDefinition("d", WomArrayType(WomStringType)) -> WomArray(WomArrayType(WomStringType), List(WomString("helloD0"), WomString("helloD1"))),
    RequiredInputDefinition("e", WomCompositeType(Map("fa" -> WomSingleFileType, "fb" -> WomIntegerType)))
      -> WomObject.withType(Map("fa" -> WomSingleFile("helloEfa"), "fb" -> WomInteger(0)), WomCompositeType(Map("fa" -> WomSingleFileType, "fb" -> WomIntegerType))),
    RequiredInputDefinition("f", WomBooleanType) -> WomBoolean.True,
    RequiredInputDefinition("g", WomBooleanType) -> WomBoolean.False,
    RequiredInputDefinition("h", WomIntegerType) -> WomInteger(0),
    RequiredInputDefinition("i", WomSingleFileType) -> WomSingleFile("ifile")
  )
  
  val localNameValues = inputs.map({ case (k, v) => k.localName -> v })
  
  val noIoFunctionSet = new IoFunctionSet {
    override def stdout(params: Seq[Try[WomValue]]): Try[WomFile] = ???
    override def glob(pattern: String): Future[Seq[String]] = ???
    override def stderr(params: Seq[Try[WomValue]]): Try[WomFile] = ???
    override def size(params: Seq[Try[WomValue]]): Future[WomFloat] = ???
    override def readFile(path: String): Future[String] = ???
    override def writeFile(path: String, content: String): Future[WomFile] = ???
  }
  
  val runtimeEnv = RuntimeEnvironment("", "", 0, 0D, 0L, 0L)
  
  def validate(tool: String, expectation: List[String]) = {
    val cltFile = better.files.File.newTemporaryFile()().write(tool)

    val clt = CwlDecoder.decodeTopLevelCwl(cltFile, None).value.unsafeRunSync() match {
      case Left(errors) => fail("Cannot parse command line tool " + errors.toList.mkString(", "))
      case Right(cwlFile) => cwlFile.select[CommandLineTool].get
    }
    
    clt
      .buildCommandTemplate(inputs)
      .valueOr(errors => fail(errors.toList.mkString(", ")))
      .flatTraverse[ErrorOr, String](_.instantiate(localNameValues, noIoFunctionSet, identity[WomValue], runtimeEnv).map(_.map(_.commandString)))
      .valueOr(errors => fail(errors.toList.mkString(", "))) shouldBe expectation
    
    cltFile.delete(true)
  }
  
  it should "filter input arguments when not bound to command line" in {
    val tool = """
      |class: CommandLineTool
      |cwlVersion: v1.0
      |inputs:
      |  - id: a
      |    type: string
      |    inputBinding: { position: 2 }
      |
      |  - id: b
      |    type: string
      |
      |baseCommand: echo
      |outputs: []
    """.stripMargin

    validate(tool, List("echo", "helloA"))
  }

  it should "order input arguments that are bound to command line" in {
    val tool = """
                 |class: CommandLineTool
                 |cwlVersion: v1.0
                 |inputs:
                 |  - id: b
                 |    type: string
                 |    inputBinding: { position: 2 }
                 |    
                 |  - id: c
                 |    type: string
                 |    inputBinding: { position: 3 }
                 |
                 |  - id: a
                 |    type: string
                 |    inputBinding: { position: 1 }
                 |
                 |baseCommand: echo
                 |outputs: []
               """.stripMargin

    validate(tool, List("echo", "helloA", "helloB", "helloC"))
  }

  it should "include arguments and sort them" in {
    val tool = """
                 |class: CommandLineTool
                 |cwlVersion: v1.0
                 |inputs:
                 |  - id: b
                 |    type: string
                 |    inputBinding: { position: 4 }
                 |    
                 |  - id: c
                 |    type: string
                 |    inputBinding: { position: 5 }
                 |
                 |  - id: a
                 |    type: string
                 |    inputBinding: { position: 1 }
                 |
                 |baseCommand: echo
                 |arguments:
                 |  - arg0
                 |  - valueFrom: "arg2"
                 |    position: 3
                 |  - valueFrom: "arg1"
                 |    position: 2
                 |outputs: []
               """.stripMargin

    validate(tool, List("echo", "arg0", "helloA", "arg1", "arg2", "helloB", "helloC"))
  }

  it should "handle array types" in {
    val tool = """
                 |class: CommandLineTool
                 |cwlVersion: v1.0
                 |inputs:
                 |  - id: d
                 |    type:
                 |      type: array
                 |      items: string
                 |      inputBinding: { prefix: "-YYY" }
                 |    inputBinding: { position: 3, prefix: "-XXX" }
                 |    
                 |  - id: a
                 |    type: string
                 |    inputBinding: { position: 1 }
                 |
                 |baseCommand: echo
                 |arguments:
                 |  - valueFrom: "arg2"
                 |    position: 2
                 |outputs: []
               """.stripMargin

    validate(tool, List("echo", "helloA", "arg2", "-XXX", "-YYY", "helloD0", "-YYY", "helloD1"))
  }

  it should "handle record schema" in {
    val tool = """
                 |class: CommandLineTool
                 |cwlVersion: v1.0
                 |inputs:
                 |  e:
                 |    type:
                 |      name: e
                 |      type: record
                 |      fields:
                 |      - name: fa
                 |        type: File
                 |        inputBinding:
                 |          position: 2
                 |          prefix: "--prefixFa"
                 |      - name: fb
                 |        type: int
                 |        inputBinding:
                 |          position: 4
                 |          prefix: "--prefixFb"
                 |    inputBinding:      
                 |      prefix: "--prefixRecord"
                 |  a:
                 |    type: string
                 |    inputBinding: { position: 1 }
                 |
                 |baseCommand: echo
                 |arguments:
                 |  - valueFrom: "arg2"
                 |    position: 3
                 |outputs: []
               """.stripMargin

    validate(tool, List("echo", "--prefixRecord", "--prefixFa", "helloEfa", "--prefixFb", "0", "helloA", "arg2"))
  }

  it should "handle prefixes properly for primitive types" in {
    val tool = """
                 |class: CommandLineTool
                 |cwlVersion: v1.0
                 |inputs:
                 |  a:
                 |    type: string
                 |    inputBinding: { position: 1, prefix: --aprefix }
                 |  f:
                 |    type: boolean
                 |    inputBinding: { position: 2, prefix: --show }
                 |  g:
                 |    type: boolean
                 |    inputBinding: { position: 3, prefix: --noshow }
                 |  h:
                 |    type: string
                 |    inputBinding: { position: 4, prefix: --hprefix }
                 |  i:
                 |    type: File
                 |    inputBinding: { position: 5, prefix: --iprefix }
                 |
                 |baseCommand: echo
                 |outputs: []
               """.stripMargin

    validate(tool, List("echo", "--aprefix", "helloA", "--show", "--hprefix", "0", "--iprefix", "ifile"))
  }

  it should "handle arrays with item separator" in {
    val tool = """
                 |class: CommandLineTool
                 |cwlVersion: v1.0
                 |inputs:
                 |  d:
                 |    type:
                 |      type: array
                 |      items: string
                 |      inputBinding: { position: 2, prefix: --item }
                 |    inputBinding: { position: 2, prefix: --array, itemSeparator: ", " }
                 |
                 |baseCommand: echo
                 |outputs: []
               """.stripMargin

    validate(tool, List("echo", "--array", "helloD0, helloD1"))
  }

  it should "handle arrays with valueFrom" in {
    val tool = """
                 |class: CommandLineTool
                 |cwlVersion: v1.0
                 |inputs:
                 |  d:
                 |    type:
                 |      type: array
                 |      items: string
                 |      inputBinding: { position: 2, prefix: --item }
                 |    inputBinding: { position: 2, prefix: --array, valueFrom: "hello" }
                 |
                 |baseCommand: echo
                 |outputs: []
               """.stripMargin

    validate(tool, List("echo", "--array", "hello"))
  }

  it should "handle the separate field" in {
    val tool = """
                 |class: CommandLineTool
                 |cwlVersion: v1.0
                 |inputs:
                 |  a:
                 |    type: string
                 |    inputBinding: { position: 0, prefix: --prefix, separate: false }
                 |
                 |baseCommand: echo
                 |outputs: []
               """.stripMargin

    validate(tool, List("echo", "--prefixhelloA"))
  }
}
