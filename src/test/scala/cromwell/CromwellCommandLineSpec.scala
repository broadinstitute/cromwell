package cromwell

import better.files._
import cromwell.core.path.PathImplicits._
import cromwell.util.SampleWdl
import cromwell.util.SampleWdl.{FileClobber, FilePassingWorkflow, ThreeStep}
import org.scalatest.{FlatSpec, Matchers}

import scala.util.Try

class CromwellCommandLineSpec extends FlatSpec with Matchers {
  import CromwellCommandLineSpec._

  behavior of "CromwellCommandLine"

  it should "UsageAndExit with no arguments" in {
    CromwellCommandLine(List.empty[String]) shouldBe UsageAndExit
  }

  it should "RunServer when specified" in {
    CromwellCommandLine(List("server")) shouldBe RunServer
  }

  it should "UsageAndExit when supplying an argument to server" in {
    CromwellCommandLine(List("server", "foo")) shouldBe UsageAndExit
  }

  it should "UsageAndExit with no arguments to run" in {
    CromwellCommandLine(List("run")) shouldBe UsageAndExit
  }

  it should "fail with too many arguments to run" in {
    CromwellCommandLine(List("run", "bork", "bork", "bork", "bork", "bork", "blerg"))
  }

  it should "VersionAndExit when the `-version` flag is passed" in {
    CromwellCommandLine(List("-version")) shouldBe VersionAndExit
  }

  it should "RunSingle when supplying wdl and inputs" in {
    CromwellCommandLine(List("run", ThreeStepWithoutOptions.wdl, ThreeStepWithoutOptions.inputs)) shouldBe a [RunSingle]
  }

  it should "RunSingle with default inputs when only supplying wdl" in {
    val ccl = CromwellCommandLine(List("run", ThreeStepWithoutOptions.wdl)).asInstanceOf[RunSingle]
     ccl.sourceFiles.inputsJson shouldBe ThreeStepInputs
  }

  it should "RunSingle with defaults if you use dashes" in {
    val ccl = CromwellCommandLine(List("run", ThreeStepWithoutOptions.wdl, "-", "-", "-")).asInstanceOf[RunSingle]
    ccl.sourceFiles.inputsJson shouldBe "{}"
    ccl.sourceFiles.workflowOptionsJson shouldBe "{}"
  }

  it should "RunSingle with options, if passed in" in {
    val threeStep = WdlAndInputs(ThreeStep, optionsJson = """{ foobar bad json! }""")
    val ccl = CromwellCommandLine(List("run", threeStep.wdl, threeStep.inputs, threeStep.options)).asInstanceOf[RunSingle]
    ccl.sourceFiles.workflowOptionsJson shouldBe threeStep.optionsJson
  }

  it should "fail if inputs path does not exist" in {
    val ccl = Try(CromwellCommandLine(List("run", ThreeStepWithoutOptions.wdl, "/some/path/that/doesnt/exit")))
    ccl.isFailure shouldBe true
    ccl.failed.get.getMessage should include("Inputs does not exist")
  }

  it should "fail if inputs path is not writeable" in {
    val threeStep = WdlAndInputs(ThreeStep)
    threeStep.inputsFile setPermissions Set.empty
    val ccl = Try(CromwellCommandLine(List("run", threeStep.wdl, threeStep.inputs)))
    ccl.isFailure shouldBe true
    ccl.failed.get.getMessage should include("Inputs is not readable")
  }

  it should "fail if metadata path is not writeable" in {
    val threeStep = WdlAndInputs(ThreeStep)
    threeStep.metadataFile write "foo"
    threeStep.metadataFile setPermissions Set.empty
    val ccl = Try(CromwellCommandLine(List("run", threeStep.wdl, threeStep.inputs, "-", threeStep.metadata)))
    ccl.isFailure shouldBe true
    ccl.failed.get.getMessage should include("Unable to write to metadata directory:")
  }

  it should "run if imports directory is a .zip file" in {
    val wdlDir = File.newTemporaryDirectory("wdlDirectory")

    val filePassing = File.newTemporaryFile("filePassing", ".wdl", Option(wdlDir))
    val fileClobber = File.newTemporaryFile("fileClobber", ".wdl", Option(wdlDir))
    filePassing write FilePassingWorkflow.wdlSource()
    fileClobber write FileClobber.wdlSource()

    val zippedDir = wdlDir.zip()
    val zippedPath = zippedDir.pathAsString

    val ccl = Try(CromwellCommandLine(List("run", filePassing.pathAsString, "-", "-", "-", zippedPath)))
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
    private var tempFiles = Vector.empty[File]

    lazy val wdlFile = {
      val file = File.newTemporaryFile(s"${sampleWdl.name}.", ".wdl")
      tempFiles :+= file
      file write sampleWdl.wdlSource("")
    }

    lazy val wdl = wdlFile.pathAsString

    lazy val inputsFile = {
      val file = File(wdlFile.path.swapExt(".wdl", ".inputs"))
      tempFiles :+= file
      file write sampleWdl.wdlJson
    }

    lazy val inputs = inputsFile.pathAsString

    lazy val optionsFile = {
      val file = File(wdlFile.path.swapExt(".wdl", ".options"))
      tempFiles :+= file
      file write optionsJson
    }

    lazy val options = optionsFile.pathAsString

    lazy val metadataFile = {
      val path = File(wdlFile.path.swapExt(".wdl", ".metadata.json"))
      tempFiles :+= path
      path
    }

    lazy val metadata = metadataFile.pathAsString

    def deleteTempFiles() = tempFiles.foreach(_.delete(swallowIOExceptions = true))
  }
}
