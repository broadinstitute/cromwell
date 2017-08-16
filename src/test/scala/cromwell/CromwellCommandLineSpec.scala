package cromwell

import cromwell.CommandLineParser.{CommandLineArguments, Run, Server}
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
    parser = CommandLineParser.buildParser()
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
    val optionsLast = parser.parse(Array("run", "3step.wdl", "--inputs", "3step.inputs"), CommandLineArguments()).get
    optionsLast.command shouldBe Some(Run)
    optionsLast.workflowSource.get.pathAsString shouldBe "3step.wdl"
    optionsLast.workflowInputs.get.pathAsString shouldBe "3step.inputs"

    val optionsFirst = parser.parse(Array("run", "--inputs", "3step.inputs", "3step.wdl"), CommandLineArguments()).get
    optionsFirst.command shouldBe Some(Run)
    optionsFirst.workflowSource.get.pathAsString shouldBe "3step.wdl"
    optionsFirst.workflowInputs.get.pathAsString shouldBe "3step.inputs"
  }

  it should "run single when supplying wdl and inputs and options" in {
    val optionsLast = parser.parse(Array("run", "3step.wdl", "--inputs", "3step.inputs", "--options", "3step.options"), CommandLineArguments()).get
    optionsLast.command shouldBe Some(Run)
    optionsLast.workflowSource.get.pathAsString shouldBe "3step.wdl"
    optionsLast.workflowInputs.get.pathAsString shouldBe "3step.inputs"
    optionsLast.workflowOptions.get.pathAsString shouldBe "3step.options"

    val optionsFirst = parser.parse(Array("run", "--inputs", "3step.inputs", "--options", "3step.options", "3step.wdl"), CommandLineArguments()).get
    optionsFirst.command shouldBe Some(Run)
    optionsFirst.workflowSource.get.pathAsString shouldBe "3step.wdl"
    optionsFirst.workflowInputs.get.pathAsString shouldBe "3step.inputs"
    optionsFirst.workflowOptions.get.pathAsString shouldBe "3step.options"
  }

  it should "fail if input files do not exist" in {
    val parsedArgs = parser.parse(Array("run", "3step.wdl", "--inputs", "3step.inputs", "--options", "3step.options"), CommandLineArguments()).get
    val validation = Try(CromwellEntryPoint.validateRunArguments(parsedArgs))

    validation.isFailure shouldBe true
    validation.failed.get.getMessage should include("Workflow source does not exist")
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
