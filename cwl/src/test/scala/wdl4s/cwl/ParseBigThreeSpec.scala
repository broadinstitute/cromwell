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

  decodeCwl(firstTool).isValid shouldBe true
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
- run: cwl/src/test/resources/arguments.cwl
  in:
  - source: untar/example_out
    id: compile/src
  out:
  - compile/classfile
  id: compile
- run: cwl/src/test/resources/tar-param.cwl
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
      .isValid should be (true)
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

    decodeCwl(envCwl).isValid should be (true)
  }
}
