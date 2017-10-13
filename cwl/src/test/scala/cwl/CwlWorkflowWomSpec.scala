package cwl

import cats.syntax.either._
import cwl.CwlDecoder._
import cwl.ExpressionEvaluator._
import eu.timepit.refined._
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import shapeless._
import wom.WomMatchers._
import wom.callable.Callable.RequiredInputDefinition
import wom.callable.{Callable, TaskDefinition, WorkflowDefinition}
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.types.{WdlFileType, WdlStringType, WdlType}

class CwlWorkflowWomSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {
  import TestSetup._

  "munging the runnable id" should "remove the filename" in {
    val id = "file:///home/dan/common-workflow-language/v1.0/examples/tar-param.cwl#example_out"
    val out = RunOutputsToTypeMap.mungeId(id)

    out shouldBe "example_out"
  }

  "munging runnable output id " should "be able to skip the path args" in {
    val id = "file:///home/dan/common-workflow-language/v1.0/examples/tar-param.cwl#ps/0b4ba500-5584-4fed-a831-9fa6f914ad3f/ps-stdOut"
    val out = RunOutputsToTypeMap.mungeId(id)

    out shouldBe "ps-stdOut"
  }


  "A Cwl object for 1st-tool" should "convert to WOM" in {
    def validateWom(callable: Callable) = callable match {
      case taskDefinition: TaskDefinition =>
        taskDefinition.inputs shouldBe List(RequiredInputDefinition(s"message", WdlStringType))
        ()

      case _ => fail("not a task definition")
    }

    (for {
      clt <- decodeAllCwl(rootPath/"1st-tool.cwl").
              map(_.select[CommandLineTool].get).
              value.
              unsafeRunSync
    } yield validateWom(clt.taskDefinition)).leftMap(e => throw new RuntimeException(s"error! $e"))
  }

  "Cwl for 1st workflow" should "convert to WOM" in {
    (for {
      wf <- decodeAllCwl(rootPath/"1st-workflow.cwl").
              value.
              unsafeRunSync.
              map(_.select[Workflow].get)

      womDefinition <- wf.womDefinition
    } yield validateWom(womDefinition)).leftMap(e => throw new RuntimeException(s"error! $e"))

    def shouldBeRequiredGraphInputNode(node: GraphNode, localName: String, wdlType: WdlType): Unit = {
      node.isInstanceOf[RequiredGraphInputNode] shouldBe true
      val requiredGraphInputNode = node.asInstanceOf[RequiredGraphInputNode]
      requiredGraphInputNode.localName shouldBe localName
      requiredGraphInputNode.womType shouldBe wdlType
      ()
    }

    def validateWom(callable: Callable) = {
      callable match {
        case wf: WorkflowDefinition =>
          val nodes = wf.innerGraph.nodes

          nodes collect {
            case gin: GraphInputNode => gin.localName
          } should be(Set("ex", "inp"))

          nodes collect {
            case cn: CallNode => cn.localName
          } should be(Set("compile", "untar"))

          val untarUpstream = nodes.collectFirst {
            case tarParam: CallNode if tarParam.localName == s"untar" => tarParam
          }.get.
            upstream
          
          untarUpstream should have size 2
          untarUpstream.collectFirst({
            case exprNode: ExpressionNode if exprNode.localName == s"file://$rootPath/1st-workflow.cwl#untar/extractfile" =>
              shouldBeRequiredGraphInputNode(exprNode.inputPorts.head.upstream.graphNode, "ex", WdlStringType)
          }).getOrElse(fail("Can't find expression node for ex"))
          
          untarUpstream.collectFirst({
            case exprNode: ExpressionNode if exprNode.localName == s"file://$rootPath/1st-workflow.cwl#untar/tarfile" =>
              shouldBeRequiredGraphInputNode(exprNode.inputPorts.head.upstream.graphNode, "inp", WdlFileType)
          }).getOrElse(fail("Can't find expression node for inp"))

          val compileUpstreamExpressionPort = nodes.collectFirst {
            case compile: CallNode if compile.localName == s"compile" => compile
          }.get.inputPorts.map(_.upstream).head

          compileUpstreamExpressionPort.name shouldBe s"file://$rootPath/1st-workflow.cwl#compile/src"
          compileUpstreamExpressionPort.graphNode.asInstanceOf[ExpressionNode].inputPorts.head.upstream.name shouldBe s"example_out"

          nodes.collect {
            case c: PortBasedGraphOutputNode => c
          }.map(_.localName) shouldBe Set(s"file://$rootPath/1st-workflow.cwl#classout")
        case wth: Any => fail(s"Parsed unexpected Callable: $wth")
      }
    }
  }

  behavior of "A decoded CWL 3step"

  private val stringOrExpressionTests = Table(
    ("index", "result"),
    (0, Coproduct[StringOrExpression]("grep")),
    (1, Coproduct[StringOrExpression](Coproduct[Expression](refineMV[MatchesECMAScript]("$(inputs.pattern)")))),
    (2, Coproduct[StringOrExpression](Coproduct[Expression](refineMV[MatchesECMAFunction]("$" + "{return inputs.file}")))),
    (3, Coproduct[StringOrExpression]("|")),
    (4, Coproduct[StringOrExpression]("wc")),
    (5, Coproduct[StringOrExpression]("-l"))
  )

  private def getTestName(stringOrExpression: StringOrExpression): String = {
    stringOrExpression.fold(StringOrExpressionToTestName)
  }

  private lazy val commandLineTool: CommandLineTool = {
    val wf = decodeAllCwl(rootPath / "three_step.cwl").map {
      _.select[Workflow].get
    }.value.unsafeRunSync.fold(error => throw new RuntimeException(s"broken parse! msg was $error"), identity)

    // The second step (aka 1) step should be cgrep
    val run: WorkflowStep.Run = wf.steps.apply(1).run
    val commandLineTool: CommandLineTool = run.select[CommandLineTool].getOrElse(fail(s"$run wasn't a CommandLineTool"))

    commandLineTool.id.get should include("cgrep")
    commandLineTool
  }

  forAll(stringOrExpressionTests) { (index, expected) =>
    it should s"correctly identify the ${getTestName(expected)}" in {
      val argument: CommandLineTool.Argument = commandLineTool.arguments.get.apply(index)
      val commandLineBinding: CommandLineBinding = argument.select[CommandLineBinding]
        .getOrElse(fail(s"$argument wasn't a CommandLineBinding"))
      val stringOrExpression: StringOrExpression = commandLineBinding.valueFrom
        .getOrElse(fail(s"valueFrom missing in $commandLineBinding"))
      stringOrExpression should be(expected)
    }
  }

  "A CwlNamespace for 3step" should "provide conversion to WOM" in {

    val wf = decodeAllCwl(rootPath/"three_step.cwl").map {
      _.select[Workflow].get
    }.value.unsafeRunSync.fold(error => throw new RuntimeException(s"broken parse! msg was $error"), identity)

    val wfd = wf.womDefinition match {
      case Right(wf: WorkflowDefinition) => wf
      case Left(o) => fail(s"Workflow definition was not produced correctly: ${o.toList.mkString(", ")}")
      case Right(callable) => fail(s"produced $callable when a Workflow Definition was expected!")
    }

    val nodes = wfd.innerGraph.nodes

    val graphInputNodes = nodes collect { case gin: GraphInputNode => gin }
    graphInputNodes should have size 1
    val patternInputNode = graphInputNodes.head
    patternInputNode.localName should be("pattern")

    nodes collect { case gon: GraphOutputNode => gon.localName } should be(Set(
      "file:///Users/danb/wdl4s/r.cwl#cgrep-count",
      "file:///Users/danb/wdl4s/r.cwl#wc-count"
    ))

    nodes collect { case cn: CallNode => cn.localName } should be(Set("ps", "cgrep", "wc"))

    val ps = nodes.collectFirst({ case ps: CallNode if ps.localName == "ps" => ps }).get
    val cgrep = nodes.collectFirst({ case cgrep: CallNode if cgrep.localName == "cgrep" => cgrep }).get
    val cgrepFileExpression = nodes.collectFirst({ case cgrepInput: ExpressionNode if cgrepInput.localName == "file:///Users/danb/wdl4s/r.cwl#cgrep/file" => cgrepInput }).get
    val cgrepPatternExpression = nodes.collectFirst({ case cgrepInput: ExpressionNode if cgrepInput.localName == "file:///Users/danb/wdl4s/r.cwl#cgrep/pattern" => cgrepInput }).get
    val wc = nodes.collectFirst({ case wc: CallNode if wc.localName == "wc" => wc }).get
    val wcFileExpression = nodes.collectFirst({ case wcInput: ExpressionNode if wcInput.localName == "file:///Users/danb/wdl4s/r.cwl#wc/file" => wcInput }).get

    ps.upstream shouldBe empty

    cgrep.upstream should contain theSameElementsAs Set(cgrepFileExpression, cgrepPatternExpression)
    wc.upstream should contain theSameElementsAs Set(wcFileExpression)

    // Check that expressions input ports point to the right output port
    cgrepPatternExpression.inputPorts.head.upstream should be theSameInstanceAs patternInputNode.singleOutputPort
    cgrepFileExpression.inputPorts.head.upstream should be theSameInstanceAs ps.outputPorts.head
    wcFileExpression.inputPorts.head.upstream should be theSameInstanceAs ps.outputPorts.head
    
    // Check that the inputDefinitionMappings are correct
    ps.inputDefinitionMappings shouldBe empty
    cgrep.inputDefinitionMappings should have size 2
    
    val cgrepFileInputDef = cgrep.callable.inputs.find(_.name == "file").get
    cgrep.inputDefinitionMappings(cgrepFileInputDef).select[OutputPort].get should be theSameInstanceAs cgrepFileExpression.singleExpressionOutputPort
    
    val cgrepPatternInputDef = cgrep.callable.inputs.find(_.name == "pattern").get
    cgrep.inputDefinitionMappings(cgrepPatternInputDef).select[OutputPort].get should be theSameInstanceAs cgrepPatternExpression.singleExpressionOutputPort
  }

}
object ExpressionTestValue extends Poly1 {
  implicit def script = at[ECMAScriptExpression] {_.value}
  implicit def function = at[ECMAScriptFunction] {_.value}
}

object StringOrExpressionToTestName extends Poly1 {
  implicit def caseECMAScript: Case.Aux[Expression, String] = {
    at[Expression] { _.fold(ExpressionTestValue) }
  }

  implicit def caseString: Case.Aux[String, String] = {
    at[String] { string => s"string $string" }
  }

}

