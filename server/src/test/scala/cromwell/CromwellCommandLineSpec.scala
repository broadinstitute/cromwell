package cromwell

import cromwell.CromwellApp.{Run, Server}
import cromwell.CromwellCommandLineSpec.WdlAndInputs
import cromwell.core.path.{DefaultPathBuilder, Path}
import cromwell.util.SampleWdl
import cromwell.util.SampleWdl.{FileClobber, FilePassingWorkflow, ThreeStep}
import org.scalatest.{BeforeAndAfter, FlatSpec, Matchers}

import scala.util.Try

class CromwellCommandLineSpec extends FlatSpec with Matchers with BeforeAndAfter {

  var parser: scopt.OptionParser[CommandLineArguments] = _

  behavior of "CromwellCommandLine"

  before {
    parser = CromwellApp.buildParser()
  }

  it should "fail to parse with no arguments" in {
    parser.parse(Array.empty[String], CommandLineArguments()).get.command shouldBe None
  }

  it should "run server when specified" in {
    parser.parse(Array("server"), CommandLineArguments()).get.command shouldBe Some(Server)
  }

  it should "fail to parse when supplying an argument to server" in {
    parser.parse(Array("server", "foo"), CommandLineArguments()) shouldBe None
  }

  it should "fail to parse with no arguments to run" in {
    parser.parse(Array("run"), CommandLineArguments()) shouldBe None
  }

  it should "fail to parse with too many arguments to run" in {
    parser.parse(Array("run", "forrest", "run"), CommandLineArguments()) shouldBe None
  }

  // --version exits the JVM which is not great in a test suite.  Haven't figure out a way to test this yet.
  //  it should "handle version output when the `-version` flag is passed" in {
  //    // I don't see a way to see that --version is printing just the version, but this at least confirms a `None`
  //    // output that should generate a usage and version.
  //    parser.parse(Array("--version"), CommandLineArguments()) shouldBe None
  //  }

  it should "run single when supplying wdl and inputs" in {
    val threeStep = WdlAndInputs(ThreeStep)
    val optionsLast = parser.parse(Array("run", threeStep.wdl, "--inputs", threeStep.inputs), CommandLineArguments()).get
    optionsLast.command shouldBe Some(Run)
    optionsLast.workflowSource.get shouldBe threeStep.wdl
    optionsLast.workflowInputs.get.pathAsString shouldBe threeStep.inputs

    val optionsFirst = parser.parse(Array("run", "--inputs", threeStep.inputs, threeStep.wdl), CommandLineArguments()).get
    optionsFirst.command shouldBe Some(Run)
    optionsFirst.workflowSource.get shouldBe threeStep.wdl
    optionsFirst.workflowInputs.get.pathAsString shouldBe threeStep.inputs

    val validation = Try(CromwellEntryPoint.validateRunArguments(optionsFirst))
    validation.isSuccess shouldBe true
    validation.get.workflowSource shouldBe None
    validation.get.workflowUrl shouldBe Some(threeStep.wdl)
  }

  it should "run single when supplying workflow using url" in {
    val url = "https://path_to_url"
    val command = parser.parse(Array("run", url), CommandLineArguments()).get

    command.command shouldBe Some(Run)
    command.workflowSource.get shouldBe url

    val validation = Try(CromwellEntryPoint.validateRunArguments(command))
    validation.isSuccess shouldBe true
    validation.get.workflowSource shouldBe None
    validation.get.workflowUrl shouldBe Option(url)
  }

  it should "run single when supplying workflow using url with inputs" in {
    val threeStep = WdlAndInputs(ThreeStep)
    val url = "https://path_to_url"
    val command = parser.parse(Array("run", url, "--inputs", threeStep.inputs), CommandLineArguments()).get

    command.command shouldBe Some(Run)
    command.workflowSource.get shouldBe url
    command.workflowInputs.get.pathAsString shouldBe threeStep.inputs

    val validation = Try(CromwellEntryPoint.validateRunArguments(command))
    validation.isSuccess shouldBe true
    validation.get.workflowSource shouldBe None
    validation.get.workflowUrl shouldBe Option(url)
  }

  it should "run single when supplying cwl workflow" in {
    val source = raw"""{"cwlVersion":"v1.0","class":"CommandLineTool","requirements":[{"class":"InlineJavascriptRequirement"}],"hints":[{"dockerPull":"debian:stretch-slim","class":"DockerRequirement"}],"inputs":[],"baseCommand":["touch","z","y","x","w","c","b","a"],"outputs":[{"type":"string","outputBinding":{"glob":"?","outputEval":"$${ return self.sort(function(a,b) { return a.location > b.location ? 1 : (a.location < b.location ? -1 : 0) }).map(f => f.basename).join(\" \") }\n""""
    val command = parser.parse(Array("run", "server/src/test/resources/cwl_glob_sort.cwl", "--type", "CWL"), CommandLineArguments()).get

    command.command shouldBe Some(Run)
    command.workflowSource.get shouldBe "server/src/test/resources/cwl_glob_sort.cwl"

    val validation = Try(CromwellEntryPoint.validateRunArguments(command))
    validation.isSuccess shouldBe true
    validation.get.workflowSource.get should include(source)
    validation.get.workflowUrl shouldBe None
  }

  it should "run single when supplying cwl workflow using url" in {
    val url = "https://path_to_url"
    val command = parser.parse(Array("run", url, "--type", "CWL"), CommandLineArguments()).get

    command.command shouldBe Some(Run)
    command.workflowSource.get shouldBe url

    val validation = Try(CromwellEntryPoint.validateRunArguments(command))
    validation.isSuccess shouldBe true
    validation.get.workflowSource shouldBe None
    validation.get.workflowUrl shouldBe Option(url)
  }

  it should "run single when supplying wdl and inputs and options" in {
    val optionsLast = parser.parse(Array("run", "3step.wdl", "--inputs", "3step.inputs", "--options", "3step.options"), CommandLineArguments()).get
    optionsLast.command shouldBe Some(Run)
    optionsLast.workflowSource.get shouldBe "3step.wdl"
    optionsLast.workflowInputs.get.pathAsString shouldBe "3step.inputs"
    optionsLast.workflowOptions.get.pathAsString shouldBe "3step.options"

    val optionsFirst = parser.parse(Array("run", "--inputs", "3step.inputs", "--options", "3step.options", "3step.wdl"), CommandLineArguments()).get
    optionsFirst.command shouldBe Some(Run)
    optionsFirst.workflowSource.get shouldBe "3step.wdl"
    optionsFirst.workflowInputs.get.pathAsString shouldBe "3step.inputs"
    optionsFirst.workflowOptions.get.pathAsString shouldBe "3step.options"
  }

  it should "fail if workflow url is invalid" in {
    val command = parser.parse(Array("run", "htpps://url_with_invalid_protocol"), CommandLineArguments()).get
    val validation = Try(CromwellEntryPoint.validateRunArguments(command))

    validation.isFailure shouldBe true
    validation.failed.get.getMessage should include("Error while validating workflow url")
  }

  it should "fail if input files do not exist" in {
    val parsedArgs = parser.parse(Array("run", "xyzshouldnotexist.wdl", "--inputs", "xyzshouldnotexist.inputs", "--options", "xyzshouldnotexist.options"), CommandLineArguments()).get
    val validation = Try(CromwellEntryPoint.validateRunArguments(parsedArgs))

    validation.isFailure shouldBe true
    validation.failed.get.getMessage should include("Workflow source path does not exist")
    validation.failed.get.getMessage should include("Workflow inputs does not exist")
    validation.failed.get.getMessage should include("Workflow options does not exist")
  }

  it should "fail if inputs path is not readable" in {
    val threeStep = WdlAndInputs(ThreeStep)
    val parsedArgs = parser.parse(Array("run", threeStep.wdl, "--inputs", threeStep.inputs), CommandLineArguments()).get
    threeStep.inputsFile setPermissions Set.empty
    val ccl = Try(CromwellEntryPoint.validateRunArguments(parsedArgs))
    ccl.isFailure shouldBe true
    ccl.failed.get.getMessage should include("Workflow inputs is not readable")
  }

  it should "fail if metadata output path is not writeable" in {
    val threeStep = WdlAndInputs(ThreeStep)
    val parsedArgs = parser.parse(Array("run", threeStep.wdl, "--inputs", threeStep.inputs, "--metadata-output", threeStep.metadata), CommandLineArguments()).get
    threeStep.metadataFile write "foo"
    threeStep.metadataFile setPermissions Set.empty
    val ccl = Try(CromwellEntryPoint.validateRunArguments(parsedArgs))
    ccl.isFailure shouldBe true
    ccl.failed.get.getMessage should include("Unable to write to metadata directory:")
  }

  it should "run if the imports path is a .zip file" in {
    val wdlDir = DefaultPathBuilder.createTempDirectory("wdlDirectory")

    val filePassing = DefaultPathBuilder.createTempFile("filePassing", ".wdl", Option(wdlDir))
    val fileClobber = DefaultPathBuilder.createTempFile("fileClobber", ".wdl", Option(wdlDir))
    filePassing write FilePassingWorkflow.workflowSource()
    fileClobber write FileClobber.workflowSource()

    val zippedDir = wdlDir.zip()
    val zippedPath = zippedDir.pathAsString

    val parsedArgs = parser.parse(Array("run", filePassing.pathAsString, "--imports", zippedPath), CommandLineArguments()).get
    val ccl = Try(CromwellEntryPoint.validateRunArguments(parsedArgs))
    ccl.isFailure shouldBe false

    zippedDir.delete(swallowIOExceptions = true)
  }
}

object CromwellCommandLineSpec {
  val ThreeStepWithoutOptions = WdlAndInputs(ThreeStep)
  val ThreeStepInputs = ThreeStepWithoutOptions.inputsFile.contentAsString

  /**
   * Create a temporary wdl file and inputs for the sampleWdl.
   * When the various properties are lazily accessed, they are also registered for deletion after the suite completes.
   */
  case class WdlAndInputs(sampleWdl: SampleWdl, optionsJson: String = "{}") {
    // Track all the temporary files we create, and delete them after the test.
    private var tempFiles = Vector.empty[Path]

    lazy val wdlFile = {
      val file = DefaultPathBuilder.createTempFile(s"${sampleWdl.name}.", ".wdl")
      tempFiles :+= file
      file write sampleWdl.workflowSource()
    }

    lazy val wdl = wdlFile.pathAsString

    lazy val inputsFile = {
      val file = wdlFile.swapExt("wdl", "inputs")
      tempFiles :+= file
      file write sampleWdl.workflowJson
    }

    lazy val inputs = inputsFile.pathAsString

    lazy val optionsFile = {
      val file = wdlFile.swapExt("wdl", "options")
      tempFiles :+= file
      file write optionsJson
    }

    lazy val options = optionsFile.pathAsString

    lazy val metadataFile = {
      val path = wdlFile.swapExt("wdl", "metadata.json")
      tempFiles :+= path
      path
    }

    lazy val metadata = metadataFile.pathAsString

    def deleteTempFiles() = tempFiles.foreach(_.delete(swallowIOExceptions = true))
  }
}
