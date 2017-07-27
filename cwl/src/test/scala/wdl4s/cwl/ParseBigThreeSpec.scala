package wdl4s.cwl


import org.scalatest._
import CwlCodecs._


class ParseBigThreeSpec extends FlatSpec with Matchers {
  val namespace = "cwl"

  it should "parse 1st tool" in {
  val firstTool = """
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

  decodeCwl(firstTool).isRight shouldBe true
  }

  it should "parse first workflow" in {
    val firstWorkflow = """
cwlVersion: v1.0
class: Workflow
inputs:
- type: string
  id: ex
- type: File
  id: inp
outputs:
- type: File
  outputSource: compile/classfile
  id: classout
steps:
- run: arguments.cwl
  in:
  - source: untar/example_out
    id: compile/src
  out:
  - compile/classfile
  id: compile
- run: tar-param.cwl
  in:
  - source: ex
    id: extractfile
  - source: inp
    id: tarfile
  out:
  - example_out
  id: untar
"""
    decodeCwl(firstWorkflow)
      .isRight should be (true)
  }

  it should "parse env cwl" in {
    val envCwl = """
cwlVersion: v1.0
class: CommandLineTool
baseCommand: env
requirements:
- envDef:
  - envValue: $(inputs.message)
    envName: HELLO
  class: EnvVarRequirement
inputs:
- type: string
  id: file:///Users/danb/common-workflow-language/v1.0/examples/env.cwl#message
outputs: []
id: file:///Users/danb/common-workflow-language/v1.0/examples/env.cwl
name: file:///Users/danb/common-workflow-language/v1.0/examples/env.cwl
"""

    val output = decodeCwl(envCwl)
    println(s"output was $output")
    output.isRight should be (true)
  }
}
