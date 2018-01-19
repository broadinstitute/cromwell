import centaur.cwl.PAPIPreprocessor
import org.scalatest.{FlatSpec, Matchers}

class PAPIPreprocessorSpec extends FlatSpec with Matchers {
  behavior of "PAPIPreProcessor"

  def validate(result: String, expectation: String) = {
    val parsedResult = io.circe.yaml.parser.parse(result).right.get
    val parsedExpectation = io.circe.yaml.parser.parse(expectation).right.get

    // This is an actual Json comparison from circe
    parsedResult shouldBe parsedExpectation
  }

  it should "prefix files" in {
    validate(
      PAPIPreprocessor.preProcessInput(
        """
          |{
          |  "input": {
          |    "file": {"location": "whale.txt", "class": "File"},
          |    "directory": {"location": "ref.fasta", "class": "Directory"}
          |  }
          |}
          |
      """.stripMargin),
      """
        |{
        |  "input": {
        |    "file": {"location": "gs://centaur-cwl-conformance/cwl-inputs/whale.txt", "class": "File"},
        |    "directory": {"location": "gs://centaur-cwl-conformance/cwl-inputs/ref.fasta", "class": "Directory"}
        |  }
        |}
        |
      """.stripMargin
    )
  }

  it should "add default docker image" in {
    validate(
      PAPIPreprocessor.preProcessWorkflow(
        """
          |class: CommandLineTool
          |cwlVersion: v1.0
          |inputs:
          |  - id: reference
          |    type: File
          |    inputBinding: { position: 2 }
          |
          |outputs:
          |  args: string[]
          |
          |baseCommand: python
          |arguments: ["bwa", "mem"]
        """.stripMargin),
      """
        |class: CommandLineTool
        |requirements:
        |  - class: DockerRequirement
        |    dockerPull: ubuntu:latest
        |cwlVersion: v1.0
        |inputs:
        |  - id: reference
        |    type: File
        |    inputBinding: { position: 2 }
        |
        |outputs:
        |  args: string[]
        |
        |baseCommand: python
        |arguments: ["bwa", "mem"]""".stripMargin
    )
  }

  it should "add default docker image in multi tool/workflow files" in {
    validate(
      PAPIPreprocessor.preProcessWorkflow(
        """
          |#!/usr/bin/env cwl-runner
          |
          |cwlVersion: v1.0
          |$graph:
          |
          |- id: echo
          |  class: CommandLineTool
          |  inputs: []
          |  outputs: []
          |  baseCommand: "echo"
          |  arguments: ["-n", "foo"]
          |
          |- id: main
          |  class: Workflow
          |  inputs: []
          |  requirements:
          |    - class: ScatterFeatureRequirement
          |  steps:
          |    step1:
          |      scatter: [echo_in1, echo_in2]
          |      scatterMethod: flat_crossproduct
          |      in: []
          |      out: []
          |
          |  outputs: []
        """.stripMargin),
      """
        |#!/usr/bin/env cwl-runner
        |
        |cwlVersion: v1.0
        |$graph:
        |
        |- id: echo
        |  class: CommandLineTool
        |  requirements:
        |  - class: DockerRequirement
        |    dockerPull: ubuntu:latest
        |  inputs: []
        |  outputs: []
        |  baseCommand: "echo"
        |  arguments: ["-n", "foo"]
        |
        |- id: main
        |  class: Workflow
        |  inputs: []
        |  requirements:
        |    - class: ScatterFeatureRequirement
        |    - class: DockerRequirement
        |      dockerPull: ubuntu:latest
        |  steps:
        |    step1:
        |      scatter: [echo_in1, echo_in2]
        |      scatterMethod: flat_crossproduct
        |      in: []
        |      out: []
        |
        |  outputs: []
        |        """.stripMargin
    )
  }

  it should "not replace existing docker hint" in {
    val workflow = """
                     |class: CommandLineTool
                     |cwlVersion: v1.0
                     |hints:
                     |- class: DockerRequirement
                     |  dockerPull: python27:slim
                     |inputs:
                     |  - id: reference
                     |    type: File
                     |    inputBinding: { position: 2 }
                     |
                     |outputs:
                     |  args: string[]
                     |
                     |baseCommand: python
                     |arguments: ["bwa", "mem"]
                   """.stripMargin
    validate(PAPIPreprocessor.preProcessWorkflow(workflow), workflow)
  }

  it should "not replace existing docker requirement" in {
    val workflow = """
                     |class: CommandLineTool
                     |cwlVersion: v1.0
                     |requirements:
                     |- class: DockerRequirement
                     |  dockerPull: python27:slim
                     |inputs:
                     |  - id: reference
                     |    type: File
                     |    inputBinding: { position: 2 }
                     |
                     |outputs:
                     |  args: string[]
                     |
                     |baseCommand: python
                     |arguments: ["bwa", "mem"]
                   """.stripMargin
    validate(PAPIPreprocessor.preProcessWorkflow(workflow), workflow)
  }
}
