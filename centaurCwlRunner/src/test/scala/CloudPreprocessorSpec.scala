import centaur.cwl.CloudPreprocessor
import com.typesafe.config.ConfigFactory
import org.scalatest.{FlatSpec, Matchers}
import wom.util.YamlUtils

class CloudPreprocessorSpec extends FlatSpec with Matchers {
  behavior of "PAPIPreProcessor"

  val pAPIPreprocessor = new CloudPreprocessor(ConfigFactory.load(), "papi.default-input-gcs-prefix")
  
  def validate(result: String, expectation: String) = {
    val parsedResult = YamlUtils.parse(result).right.get
    val parsedExpectation = YamlUtils.parse(expectation).right.get

    // This is an actual Json comparison from circe
    parsedResult shouldBe parsedExpectation
  }

  it should "prefix files and directories in inputs" in {
    validate(
      pAPIPreprocessor.preProcessInput(
        """|{
           |  "input": {
           |    "null": null,
           |    "file": {
           |      "location": "whale.txt",
           |      "class": "File",
           |      "secondaryFiles": [
           |        {
           |          "class": File,
           |          "location": "hello.txt"
           |        }
           |      ],
           |      "default": {
           |        "location": "default_whale.txt",
           |        "class": "File",
           |      }
           |    },
           |    "directory": {
           |      "location": "ref.fasta",
           |      "class": "Directory",
           |      "listing": [
           |        {
           |          "class": File,
           |          "location": "hello.txt"
           |        }
           |      ]
           |    }
           |  }
           |}
           |""".stripMargin).value.unsafeRunSync().right.get,
      """|{
         |  "input": {
         |    "null": null,
         |    "file": {
         |      "location": "gs://centaur-cwl-conformance-1f501e3/cwl-inputs/whale.txt",
         |      "class": "File",
         |      "secondaryFiles": [
         |        {
         |          "class": File,
         |          "location": "gs://centaur-cwl-conformance-1f501e3/cwl-inputs/hello.txt"
         |        }
         |      ],
         |      "default": {
         |        "location": "gs://centaur-cwl-conformance-1f501e3/cwl-inputs/default_whale.txt",
         |        "class": "File",
         |      }
         |    },
         |    "directory": {
         |      "location": "gs://centaur-cwl-conformance-1f501e3/cwl-inputs/ref.fasta",
         |      "class": "Directory",
         |      "listing": [
         |        {
         |          "class": File,
         |          "location": "gs://centaur-cwl-conformance-1f501e3/cwl-inputs/hello.txt"
         |        }
         |      ]
         |    }
         |  }
         |}
         |""".stripMargin
    )
  }

  it should "prefix files and directories in workflow" in {
    validate(
      pAPIPreprocessor.preProcessWorkflow(
        """|class: CommandLineTool
           |cwlVersion: v1.0
           |requirements:
           |  - class: DockerRequirement
           |    dockerPull: ubuntu:latest
           |inputs:
           |  - id: reference
           |    type: File
           |    inputBinding: { position: 2 }
           |    default:
           |      class: File
           |      location: args.py
           |
           |outputs:
           |  args: string[]
           |
           |baseCommand: python
           |arguments: ["bwa", "mem"]
           |""".stripMargin),
      """|class: CommandLineTool
         |cwlVersion: v1.0
         |requirements:
         |  - class: DockerRequirement
         |    dockerPull: ubuntu:latest
         |inputs:
         |  - id: reference
         |    type: File
         |    inputBinding: { position: 2 }
         |    default:
         |      class: File
         |      location: gs://centaur-cwl-conformance-1f501e3/cwl-inputs/args.py
         |
         |outputs:
         |  args: string[]
         |
         |baseCommand: python
         |arguments: ["bwa", "mem"]
         |""".stripMargin
    )
  }

  it should "add default docker image if there's no hint" in {
    validate(
      pAPIPreprocessor.preProcessWorkflow(
        """|class: CommandLineTool
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
           |""".stripMargin),
      """|class: CommandLineTool
         |hints:
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
         |arguments: ["bwa", "mem"]
         |""".stripMargin
    )
  }

  it should "append default docker image to existing hint as an array" in {
    validate(
      pAPIPreprocessor.preProcessWorkflow(
        """|class: CommandLineTool
           |hints:
           |  - class: EnvVarRequirement
           |    envDef:
           |      TEST_ENV: $(inputs.in)
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
           |""".stripMargin),
      """|class: CommandLineTool
         |hints:
         |  - class: EnvVarRequirement
         |    envDef:
         |      TEST_ENV: $(inputs.in)
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
         |arguments: ["bwa", "mem"]
         |""".stripMargin
    )
  }

  it should "append default docker image to existing hints as an object" in {
    validate(
      pAPIPreprocessor.preProcessWorkflow(
        """|class: CommandLineTool
           |cwlVersion: v1.0
           |inputs:
           |  in: string
           |outputs:
           |  out:
           |    type: File
           |    outputBinding:
           |      glob: out
           |
           |hints:
           |  EnvVarRequirement:
           |    envDef:
           |      TEST_ENV: $(inputs.in)
           |
           |baseCommand: ["/bin/bash", "-c", "echo $TEST_ENV"]
           |
           |stdout: out
           |""".stripMargin),
      """|class: CommandLineTool
         |cwlVersion: v1.0
         |inputs:
         |  in: string
         |outputs:
         |  out:
         |    type: File
         |    outputBinding:
         |      glob: out
         |
         |hints:
         |  EnvVarRequirement:
         |    envDef:
         |      TEST_ENV: $(inputs.in)
         |  DockerRequirement:
         |    dockerPull: ubuntu:latest
         |
         |baseCommand: ["/bin/bash", "-c", "echo $TEST_ENV"]
         |
         |stdout: out
         |""".stripMargin
    )
  }

  it should "not replace existing docker requirement in an object" in {
    val workflow = """|class: CommandLineTool
                      |cwlVersion: v1.0
                      |requirements:
                      |  DockerRequirement:
                      |    dockerPull: python27:slim
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
                      |""".stripMargin
    validate(pAPIPreprocessor.preProcessWorkflow(workflow), workflow)
  }

  it should "not replace existing docker hint in an object" in {
    val workflow = """|class: CommandLineTool
                      |cwlVersion: v1.0
                      |hints:
                      |  DockerRequirement:
                      |    dockerPull: python27:slim
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
                      |""".stripMargin
    validate(pAPIPreprocessor.preProcessWorkflow(workflow), workflow)
  }

  it should "not replace existing docker hint in an array" in {
    val workflow = """|class: CommandLineTool
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
                      |""".stripMargin
    validate(pAPIPreprocessor.preProcessWorkflow(workflow), workflow)
  }

  it should "not replace existing docker requirement in an array" in {
    val workflow = """|class: CommandLineTool
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
                      |""".stripMargin
    validate(pAPIPreprocessor.preProcessWorkflow(workflow), workflow)
  }

  it should "throw an exception if yaml / json can't be parse" in {
    val invalid =
      """
        |{ [invalid]: }
      """.stripMargin

    an[Exception] shouldBe thrownBy(pAPIPreprocessor.preProcessWorkflow(invalid))
  }
}
