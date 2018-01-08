package cwl

import cats.syntax.either._
import cwl.CwlDecoder._
import cwl.ExpressionEvaluator._
import cwl.command.ParentName
import eu.timepit.refined._
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import shapeless._
import wom.WomMatchers._
import wom.callable.Callable.RequiredInputDefinition
import wom.callable.{Callable, TaskDefinition, WorkflowDefinition}
import wom.graph.GraphNodePort.OutputPort
import wom.graph._
import wom.graph.expression.ExpressionNode
import wom.types.{WomSingleFileType, WomStringType, WomType}

class CwlWorkflowWomSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {
  import TestSetup._

  implicit val parentName = ParentName.empty
  
  "A Cwl object for 1st-tool" should "convert to WOM" in {
    def validateWom(callable: Callable): Unit = callable match {
      case taskDefinition: TaskDefinition =>
        taskDefinition.inputs shouldBe List(RequiredInputDefinition(s"message", WomStringType))
        ()

      case _ => fail("not a task definition")
    }

    import cats.syntax.validated._
    (for {
      clt <- decodeAllCwl(rootPath/"1st-tool.cwl").
              map(_.select[CommandLineTool].get).
              value.
              unsafeRunSync
      taskDef <- clt.buildTaskDefinition(_.validNel).toEither
    } yield validateWom(taskDef)).leftMap(e => throw new RuntimeException(s"error! $e"))
  }

  "Cwl for 1st workflow" should "convert to WOM" in {
    (for {
      wf <- decodeAllCwl(rootPath/"1st-workflow.cwl").
              value.
              unsafeRunSync.
              map(_.select[Workflow].get)

      womDefinition <- wf.womDefinition(AcceptAllRequirements)
    } yield validateWom(womDefinition)).leftMap(e => throw new RuntimeException(s"error! ${e.toList.mkString("\n")}"))

    def shouldBeRequiredGraphInputNode(node: GraphNode, localName: String, womType: WomType): Unit = {
      node.isInstanceOf[RequiredGraphInputNode] shouldBe true
      val requiredGraphInputNode = node.asInstanceOf[RequiredGraphInputNode]
      requiredGraphInputNode.localName shouldBe localName
      requiredGraphInputNode.womType shouldBe womType
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
              shouldBeRequiredGraphInputNode(exprNode.inputPorts.head.upstream.graphNode, "ex", WomStringType)
          }).getOrElse(fail("Can't find expression node for ex"))

          untarUpstream.collectFirst({
            case exprNode: ExpressionNode if exprNode.localName == s"file://$rootPath/1st-workflow.cwl#untar/tarfile" =>
              exprNode.inputPorts.map(_.upstream.graphNode).count {
                case rgin: RequiredGraphInputNode =>
                  rgin.identifier.localName == LocalName("inp") &&
                    rgin.womType == WomSingleFileType
              }  shouldBe 1
          }).getOrElse(fail("Can't find expression node for inp"))

          val compileUpstreamExpressionPort = nodes.collectFirst {
            case compile: CallNode if compile.localName == s"compile" => compile
          }.get.inputPorts.map(_.upstream).head

          compileUpstreamExpressionPort.name shouldBe s"file://$rootPath/1st-workflow.cwl#compile/src"
          compileUpstreamExpressionPort.graphNode.asInstanceOf[ExpressionNode].inputPorts.map(_.upstream.name).count(_ == "example_out") shouldBe 1

          nodes.collect {
            case c: PortBasedGraphOutputNode => c
          }.map(_.localName) shouldBe Set("classout")
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

  private def getTestName(stringOrExpression: StringOrExpression): String = stringOrExpression match {
    case StringOrExpression.String(s) => s"string $s"
    // Although these two cases look the same, they're actually calling different 'value' functions. So we
    // can't collapse this to a single "case StringOrExpression.Expression(e) => e.value":
    case StringOrExpression.ECMAScriptExpression(e) => e.value
    case StringOrExpression.ECMAScriptFunction(f) => f.value
  }

  private lazy val commandLineTool: CommandLineTool = {
    val wf = decodeAllCwl(rootPath / "three_step.cwl").map { wf =>
      wf.select[Workflow].get
    }.value.unsafeRunSync.fold(error => throw new RuntimeException(s"broken parse! msg was ${error.toList.mkString(", ")}"), identity)

    wf.id should include("three_step")

    // The second step (aka 1) step should be cgrep
    val run: WorkflowStep.Run = wf.steps.apply(1).run
    val commandLineTool: CommandLineTool = run.select[CommandLineTool].getOrElse(fail(s"$run wasn't a CommandLineTool"))

    commandLineTool.id should include("cgrep")
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
    }.value.unsafeRunSync.fold(error => throw new RuntimeException(s"broken parse: $error"), identity)

    wf.id should include("three_step")

    val wfd = wf.womDefinition(AcceptAllRequirements) match {
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
      "cgrep-count",
      "wc-count"
    ))

    nodes collect { case cn: CallNode => cn.localName } should be(Set("ps", "cgrep", "wc"))

    val ps = nodes.collectFirst({ case ps: CallNode if ps.localName == "ps" => ps }).get
    val cgrep = nodes.collectFirst({ case cgrep: CallNode if cgrep.localName == "cgrep" => cgrep }).get
    val cgrepFileExpression = nodes.collectFirst({ case cgrepInput: ExpressionNode if s"${FileStepAndId(cgrepInput.localName).stepId}/${FileStepAndId(cgrepInput.localName).id}" == "cgrep/file" => cgrepInput }).get
    val cgrepPatternExpression = nodes.collectFirst({ case cgrepInput: ExpressionNode if s"${FileStepAndId(cgrepInput.localName).stepId}/${FileStepAndId(cgrepInput.localName).id}" == "cgrep/pattern" => cgrepInput }).get
    val wc = nodes.collectFirst({ case wc: CallNode if wc.localName == "wc" => wc }).get
    val wcFileExpression = nodes.collectFirst({ case wcInput: ExpressionNode if s"${FileStepAndId(wcInput.localName).stepId}/${FileStepAndId(wcInput.localName).id}" == "wc/file" => wcInput }).get

    ps.upstream shouldBe empty

    cgrep.upstream should contain theSameElementsAs Set(cgrepFileExpression, cgrepPatternExpression)
    wc.upstream should contain theSameElementsAs Set(wcFileExpression)

    // Check that expressions input ports point to the right output port
    cgrepPatternExpression.inputPorts.head.upstream should be theSameInstanceAs patternInputNode.singleOutputPort
    cgrepFileExpression.inputPorts.map(_.upstream).count(_ eq ps.outputPorts.head) shouldBe 1
    wcFileExpression.inputPorts.map(_.upstream).count(_ eq ps.outputPorts.head) shouldBe 1

    // Check that the inputDefinitionMappings are correct
    ps.inputDefinitionMappings shouldBe empty
    cgrep.inputDefinitionMappings should have size 2

    val cgrepInputs = cgrep.inputDefinitionMappings.toMap
    val cgrepFileInputDef = cgrep.callable.inputs.find(_.name == "file").get
    cgrepInputs(cgrepFileInputDef).select[OutputPort].get should be theSameInstanceAs cgrepFileExpression.singleExpressionOutputPort

    val cgrepPatternInputDef = cgrep.callable.inputs.find(_.name == "pattern").get
    cgrepInputs(cgrepPatternInputDef).select[OutputPort].get should be theSameInstanceAs cgrepPatternExpression.singleExpressionOutputPort
  }

}
