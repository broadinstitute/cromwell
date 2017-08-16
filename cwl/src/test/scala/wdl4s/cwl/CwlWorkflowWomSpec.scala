package wdl4s.cwl

import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.wdl.types.{WdlFileType, WdlStringType}
import wdl4s.wom.callable.Callable.RequiredInputDefinition
import wdl4s.wom.callable.{TaskDefinition, WorkflowDefinition}
import wdl4s.wom.executable.Executable
import wdl4s.wom.graph._

class CwlWorkflowWomSpec extends FlatSpec with Matchers {
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
    val firstTool =
      """
cwlVersion: v1.0
class: CommandLineTool
baseCommand: echo
inputs:
- type: string
  inputBinding:
    position: 1
  id: message
outputs: []
"""

    CwlCodecs.decodeCwl(firstTool) map {
      case (clt: CommandLineTool, _) =>
        clt.womExecutable match {
          case Valid(wom) =>
            wom.entryPoint match {
              case taskDefinition: TaskDefinition =>

                taskDefinition.inputs shouldBe Set(RequiredInputDefinition("message", WdlStringType))

                taskDefinition.graph.map {
                  graph =>
                    graph.nodes.collect { case gin: GraphInputNode => gin.name } should be(Set("echo.message"))
                    graph.nodes collect { case cn: CallNode => cn.name } should be(Set("echo"))

                    graph.nodes.collectFirst { case echo: CallNode if echo.name == "echo" => echo }.get.
                      upstream shouldBe Set(RequiredGraphInputNode("echo.message", WdlStringType))
                }

              case _ => fail("not a task definition")
            }
          case i@Invalid(_) => fail(s"Invalid: $i")
        }
      case (wth: Any, _) => fail(s"parsed unexpected type $wth")
    } leftMap { error => fail(s"did not parse!  $error") }
  }

  "Cwl for 1st workflow" should "convert to WOM" in {
    val firstWorkflow =
      s"""
cwlVersion: "v1.0"
class: "Workflow"
inputs:
  - type: "string"
    id: "file:///home/dan/wdl4s/r.cwl#ex"
  - type: "File"
    id: "file:///home/dan/wdl4s/r.cwl#inp"
outputs:
  - type: "File"
    outputSource: "file:///home/dan/wdl4s/r.cwl#compile/classfile"
    id: "file:///home/dan/wdl4s/r.cwl#classout"
steps:
  - run: "cwl/src/test/resources/arguments.cwl"
    in:
      -
        source: "file:///home/dan/wdl4s/r.cwl#untar/example_out"
        id: "file:///home/dan/wdl4s/r.cwl#compile/src"
    out:
      - "file:///home/dan/wdl4s/r.cwl#compile/classfile"
    id: "file:///home/dan/wdl4s/r.cwl#compile"
  - run: "cwl/src/test/resources/tar-param.cwl"
    in:
      -
        source: "file:///home/dan/wdl4s/r.cwl#ex"
        id: "file:///home/dan/wdl4s/r.cwl#untar/extractfile"
      -
        source: "file:///home/dan/wdl4s/r.cwl#inp"
        id: "file:///home/dan/wdl4s/r.cwl#untar/tarfile"
    out:
      - "file:///home/dan/wdl4s/r.cwl#untar/example_out"
    id: "file:///home/dan/wdl4s/r.cwl#untar"
id: "file:///home/dan/wdl4s/r.cwl"
name: "file:///home/dan/wdl4s/r.cwl"

"""

    import CwlCodecs._


    decodeCwl(firstWorkflow) map {
      case (workflow: Workflow, nameToFile) =>
        workflow.womExecutable(nameToFile) match {
          case Valid(ex) => validateWom(ex)
          case Invalid(throwable) => fail(s"executable was not created $throwable")
        }
      case (wth: Any, _) => fail(s"Parsed unexpected CwlFile subtype $wth")
    } leftMap { error => fail(s"did not parse!  $error") }

    def validateWom(ex: Executable) = {
      ex match {
        case Executable(wf: WorkflowDefinition) =>
          val nodes = wf.innerGraph.nodes

          nodes collect {
            case gin: GraphInputNode => gin.name
          } should be(Set("file:///home/dan/wdl4s/r.cwl#ex", "file:///home/dan/wdl4s/r.cwl#inp"))

          nodes collect {
            case cn: CallNode => cn.name
          } should be(Set("file:///home/dan/wdl4s/r.cwl#compile", "file:///home/dan/wdl4s/r.cwl#untar"))

          nodes.collectFirst {
            case tarParam: CallNode if tarParam.name == "file:///home/dan/wdl4s/r.cwl#untar" => tarParam
          }.get.
            upstream shouldBe Set(
            RequiredGraphInputNode("file:///home/dan/wdl4s/r.cwl#ex", WdlStringType),
            RequiredGraphInputNode("file:///home/dan/wdl4s/r.cwl#inp", WdlFileType))


          nodes.collectFirst {
            case compile: CallNode if compile.name == "file:///home/dan/wdl4s/r.cwl#compile" => compile
          }.get.inputPorts.map(_.upstream).head.name shouldBe "file:///home/dan/wdl4s/r.cwl#untar/example_out"

          nodes.collect {
            case c: PortBasedGraphOutputNode => c
          }.map(_.name) shouldBe Set("file:///home/dan/wdl4s/r.cwl#classout")
        case Executable(wth: Any) => fail(s"Parsed unexpected Callable: $wth")
      }
    }
  }

  "A WdlNamespace for 3step" should "provide conversion to WOM" in {
    val threeStep =
      """
cwlVersion: v1.0
class: Workflow
inputs:
- id: file:///Users/danb/wdl4s/r.cwl#pattern
  type: string
outputs:
- id: file:///Users/danb/wdl4s/r.cwl#cgrep-count
  outputSource: file:///Users/danb/wdl4s/r.cwl#cgrep/cgrep-count
  type: int
- id: file:///Users/danb/wdl4s/r.cwl#wc-count
  outputSource: file:///Users/danb/wdl4s/r.cwl#wc/wc-count
  type: int
steps:
- id: file:///Users/danb/wdl4s/r.cwl#ps
  in: []
  out:
  - file:///Users/danb/wdl4s/r.cwl#ps/ps-stdOut
  run:
    inputs: []
    outputs:
    - id: file:///Users/danb/wdl4s/r.cwl#ps/0b4ba500-5584-4fed-a831-9fa6f914ad3f/ps-stdOut
      outputBinding:
        glob: ps-stdOut.txt
      type: File
    class: CommandLineTool
    baseCommand: ps
    stdout: ps-stdOut.txt
    id: file:///Users/danb/wdl4s/r.cwl#ps/0b4ba500-5584-4fed-a831-9fa6f914ad3f
- id: file:///Users/danb/wdl4s/r.cwl#cgrep
  in:
  - id: file:///Users/danb/wdl4s/r.cwl#cgrep/pattern
    source: file:///Users/danb/wdl4s/r.cwl#pattern
  - id: file:///Users/danb/wdl4s/r.cwl#cgrep/file
    source: file:///Users/danb/wdl4s/r.cwl#ps/ps-stdOut
  out:
  - id: file:///Users/danb/wdl4s/r.cwl#cgrep/cgrep-count
  run:
    inputs:
    - id: file:///Users/danb/wdl4s/r.cwl#cgrep/09f8bcac-a91a-49d5-afb6-2f1b1294e875/pattern
      type: string
    - id: file:///Users/danb/wdl4s/r.cwl#cgrep/09f8bcac-a91a-49d5-afb6-2f1b1294e875/file
      type: File
    outputs:
    - id: file:///Users/danb/wdl4s/r.cwl#cgrep/09f8bcac-a91a-49d5-afb6-2f1b1294e875/cgrep-count
      outputBinding:
        glob: cgrep-stdOut.txt
        loadContents: true
        outputEval: $(self[0].contents.toInt)
      type: int
    class: CommandLineTool
    requirements:
    - class: ShellCommandRequirement
    - class: InlineJavascriptRequirement
    arguments:
    - valueFrom: grep
      shellQuote: false
    - valueFrom: $(inputs.pattern).
      shellQuote: false
    - valueFrom: $(inputs.file)
      shellQuote: false
    - valueFrom: '|'
      shellQuote: false
    - valueFrom: wc
      shellQuote: false
    - valueFrom: -l
      shellQuote: false
    stdout: cgrep-stdOut.txt
    id: file:///Users/danb/wdl4s/r.cwl#cgrep/09f8bcac-a91a-49d5-afb6-2f1b1294e875
- id: file:///Users/danb/wdl4s/r.cwl#wc
  in:
  - id: file:///Users/danb/wdl4s/r.cwl#wc/file
    source: file:///Users/danb/wdl4s/r.cwl#ps/ps-stdOut
  out:
  - id: file:///Users/danb/wdl4s/r.cwl#wc/wc-count
  run:
    inputs:
    - id: file:///Users/danb/wdl4s/r.cwl#wc/45d98851-7bfe-473e-ab24-aac922553f3e/file
      type: File
    outputs:
    - id: file:///Users/danb/wdl4s/r.cwl#wc/45d98851-7bfe-473e-ab24-aac922553f3e/wc-count
      outputBinding:
        glob: wc-stdOut.txt
        loadContents: true
        outputEval: $(self[0].contents.toInt)
      type: int
    class: CommandLineTool
    requirements:
    - class: ShellCommandRequirement
    - class: InlineJavascriptRequirement
    arguments:
    - valueFrom: cat
      shellQuote: false
    - valueFrom: $(inputs.file)
      shellQuote: false
    - valueFrom: '|'
      shellQuote: false
    - valueFrom: wc
      shellQuote: false
    - valueFrom: -l
      shellQuote: false
    stdout: wc-stdOut.txt
    id: file:///Users/danb/wdl4s/r.cwl#wc/45d98851-7bfe-473e-ab24-aac922553f3e
id: file:///Users/danb/wdl4s/r.cwl
"""

    val workflow = CwlCodecs.decodeCwl(threeStep) match {
      case Valid((wf: Workflow, _)) => wf
      case o => fail(s"didn't parse a workflow! $o")
    }
    val wf = workflow.womExecutable(Map.empty) match {
      case Valid(Executable(wf: WorkflowDefinition)) => wf
      case o => fail(s"invalid executable $o")
    }

    val nodes = wf.innerGraph.nodes

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
