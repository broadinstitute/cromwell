package cromwell

import java.nio.file.Path
import better.files._
import cromwell.core.PathFactory._
import cromwell.util.SampleWdl
import cromwell.util.SampleWdl.ThreeStep
import org.scalatest.{FlatSpec, Matchers}

import scala.language.postfixOps
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
    CromwellCommandLine(List("run", "bork", "bork", "bork", "bork", "bork"))
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
    threeStep.inputsPath setPermissions Set.empty
    val ccl = Try(CromwellCommandLine(List("run", threeStep.wdl, threeStep.inputs)))
    ccl.isFailure shouldBe true
    ccl.failed.get.getMessage should include("Inputs is not readable")
  }

  it should "fail if metadata path is not writeable" in {
    val threeStep = WdlAndInputs(ThreeStep)
    threeStep.metadataPath write "foo"
    threeStep.metadataPath setPermissions Set.empty
    val ccl = Try(CromwellCommandLine(List("run", threeStep.wdl, threeStep.inputs, "-", threeStep.metadata)))
    ccl.isFailure shouldBe true
    ccl.failed.get.getMessage should include("Unable to write to metadata directory:")
  }
}

object CromwellCommandLineSpec {
  val ThreeStepWithoutOptions = WdlAndInputs(ThreeStep)
  val ThreeStepInputs = ThreeStepWithoutOptions.inputsPath.contentAsString

  /**
   * Create a temporary wdl file and inputs for the sampleWdl.
   * When the various properties are lazily accessed, they are also registered for deletion after the suite completes.
   */
  case class WdlAndInputs(sampleWdl: SampleWdl, optionsJson: String = "{}") {
    // Track all the temporary files we create, and delete them after the test.
    private var tempFiles = Vector.empty[Path]

    lazy val wdlPath: Path = {
      val path = File.newTemp(s"${sampleWdl.name}.", ".wdl").path
      tempFiles :+= path
      path write sampleWdl.wdlSource("")
      path
    }

    lazy val wdl = wdlPath.fullPath

    lazy val inputsPath = {
      val path = wdlPath.swapExt(".wdl", ".inputs")
      tempFiles :+= path
      path write sampleWdl.wdlJson
      path
    }

    lazy val inputs = inputsPath.fullPath

    lazy val optionsPath = {
      val path = wdlPath.swapExt(".wdl", ".options")
      tempFiles :+= path
      path write optionsJson
      path
    }

    lazy val options = optionsPath.fullPath

    lazy val metadataPath = {
      val path = wdlPath.swapExt(".wdl", ".metadata.json")
      tempFiles :+= path
      path.toAbsolutePath
    }

    lazy val metadata = metadataPath.fullPath

    def deleteTempFiles() = tempFiles.foreach(_.delete(ignoreIOExceptions = true))
  }
}
