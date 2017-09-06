package wdl4s.cwl

import cats.data.Validated.Valid
import cats.syntax.either._
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.wdl.types.{WdlFileType, WdlStringType}
import wdl4s.wom.callable.Callable.RequiredInputDefinition
import wdl4s.wom.callable.{TaskDefinition, WorkflowDefinition}
import wdl4s.wom.executable.Executable
import wdl4s.wom.graph._
import CwlDecoder._

class CwlWorkflowWomSpec extends FlatSpec with Matchers {
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
    val validateWom: Executable => Unit =
      _.entryPoint match {
        case taskDefinition: TaskDefinition =>
          taskDefinition.inputs shouldBe Set(RequiredInputDefinition(s"file://$rootPath/1st-tool.cwl#message", WdlStringType))

          taskDefinition.graph.map {
            graph =>
              graph.nodes.collect { case gin: GraphInputNode => gin.name } should be(Set(s"file://$rootPath/1st-tool.cwl.file://$rootPath/1st-tool.cwl#message"))
              graph.nodes collect { case cn: CallNode => cn.name } should be(Set(s"file://$rootPath/1st-tool.cwl"))
          }
          ()

        case _ => fail("not a task definition")
      }

    (for {
      clt <- decodeAllCwl(rootPath/"1st-tool.cwl").
              map(_.select[CommandLineTool].get).
              value.
              unsafeRunSync
      wom <-  clt.womExecutable.toEither
    } yield validateWom(wom)).leftMap(e => throw new RuntimeException(s"error! $e"))
  }

  "Cwl for 1st workflow" should "convert to WOM" in {
    (for {
      wf <- decodeAllCwl(rootPath/"1st-workflow.cwl").
              value.
              unsafeRunSync.
              map(_.select[Workflow].get)

      ex <- wf.womExecutable.toEither
    } yield validateWom(ex)).leftMap(e => throw new RuntimeException(s"error! $e"))

    def validateWom(ex: Executable) = {
      ex match {
        case Executable(wf: WorkflowDefinition) =>
          val nodes = wf.innerGraph.nodes

          nodes collect {
            case gin: GraphInputNode => gin.name
          } should be(Set(s"file://$rootPath/1st-workflow.cwl#ex", s"file://$rootPath/1st-workflow.cwl#inp"))

          nodes collect {
            case cn: CallNode => cn.name
          } should be(Set(s"file://$rootPath/1st-workflow.cwl#compile", s"file://$rootPath/1st-workflow.cwl#untar"))

          nodes.collectFirst {
            case tarParam: CallNode if tarParam.name == s"file://$rootPath/1st-workflow.cwl#untar" => tarParam
          }.get.
            upstream shouldBe Set(
            RequiredGraphInputNode(s"file://$rootPath/1st-workflow.cwl#ex", WdlStringType),
            RequiredGraphInputNode(s"file://$rootPath/1st-workflow.cwl#inp", WdlFileType))


          nodes.collectFirst {
            case compile: CallNode if compile.name == s"file://$rootPath/1st-workflow.cwl#compile" => compile
          }.get.inputPorts.map(_.upstream).head.name shouldBe s"file://$rootPath/1st-workflow.cwl#untar/example_out"

          nodes.collect {
            case c: PortBasedGraphOutputNode => c
          }.map(_.name) shouldBe Set(s"file://$rootPath/1st-workflow.cwl#classout")
        case Executable(wth: Any) => fail(s"Parsed unexpected Callable: $wth")
      }
    }
  }

  "A WdlNamespace for 3step" should "provide conversion to WOM" in {

    val wf = decodeAllCwl(rootPath/"three_step.cwl").map {
      _.select[Workflow].get
    }.value.unsafeRunSync.fold(error => throw new RuntimeException(s"broken parse! msg was $error"), identity)

    val wfd = wf.womExecutable match {
      case Valid(Executable(wf: WorkflowDefinition)) => wf
      case o => fail(s"invalid executable $o")
    }

    val nodes = wfd.innerGraph.nodes

    nodes collect { case gin: GraphInputNode => gin.name } should be(Set("file:///Users/danb/wdl4s/r.cwl#pattern"))

    nodes collect { case gon: GraphOutputNode => gon.name } should be(Set(
      "file:///Users/danb/wdl4s/r.cwl#cgrep-count",
      "file:///Users/danb/wdl4s/r.cwl#wc-count"
    ))

    nodes collect { case cn: CallNode => cn.name } should be(
      Set(
        "file:///Users/danb/wdl4s/r.cwl#ps",
        "file:///Users/danb/wdl4s/r.cwl#cgrep",
        "file:///Users/danb/wdl4s/r.cwl#wc"))

    val ps = nodes.collectFirst({ case ps: CallNode if ps.name == "file:///Users/danb/wdl4s/r.cwl#ps" => ps }).get
    val cgrep = nodes.collectFirst({ case cgrep: CallNode if cgrep.name == "file:///Users/danb/wdl4s/r.cwl#cgrep" => cgrep }).get
    val cgrepPatternInput = nodes.collectFirst({ case cgrepInput: GraphInputNode if cgrepInput.name == "file:///Users/danb/wdl4s/r.cwl#pattern" => cgrepInput }).get
    val wc = nodes.collectFirst({ case wc: CallNode if wc.name == "file:///Users/danb/wdl4s/r.cwl#wc" => wc }).get

    ps.upstream shouldBe empty
    cgrep.upstream.filter(_ eq ps).size shouldBe 1
    cgrep.upstream.filter(_ eq cgrepPatternInput).size shouldBe 1

    wc.upstream.filter(_ eq ps).size shouldBe 1
  }

}
