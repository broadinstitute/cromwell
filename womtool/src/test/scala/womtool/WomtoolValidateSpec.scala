package womtool

import java.nio.file.Path

import better.files.File
import cromwell.core.path.DefaultPathBuilder
import org.scalatest.{FlatSpec, Matchers}
import womtool.WomtoolMain.{SuccessfulTermination, UnsuccessfulTermination}

import scala.collection.immutable

class WomtoolValidateSpec extends FlatSpec with Matchers {

  private val presentWorkingDirectoryName = DefaultPathBuilder.get(".").toAbsolutePath.name
  val validationTestCases = File("womtool/src/test/resources/validate")
  val languageVersions = Option(validationTestCases.list).toList.flatten

  behavior of "womtool validate"

  it should "test at least one version" in {
    languageVersions.isEmpty should be(false)
  }

  // The filterNot(_.contains(".DS")) stuff prevents Mac 'Desktop Services' hidden directories from accidentally being picked up:
  languageVersions.filterNot(f => f.name.contains(".DS")) foreach { versionDirectory =>
    val versionName = versionDirectory.name

    val validTestCases = versionDirectory.path.resolve("valid")
    val invalidTestCases = versionDirectory.path.resolve("invalid")
    val ignoredTestCases = versionDirectory.path.resolve("ignored")

    // Don't bother checking that the 'ignored' directory exists:
    List(validTestCases, invalidTestCases) foreach { path =>
      it should s"be set up for testing $versionName in '${versionDirectory.relativize(path).toString}'" in {
        if (!path.toFile.exists) fail(s"Path doesn't exist: ${path.toAbsolutePath.toString}")
        if (Option(path.toFile.list).toList.flatten.isEmpty) fail(s"No test cases found in: ${path.toAbsolutePath.toString}")
        versionDirectory.list.nonEmpty shouldBe true
      }
    }

    Option(ignoredTestCases.toFile.list).toList.flatten foreach { ignoredCase =>
      it should s"run $versionName test '$ignoredCase'" ignore {}
    }

    listFilesAndFilterDSFile(validTestCases) foreach { validCase =>
      val inputsFile = ifExists(validTestCases.resolve(validCase).resolve(validCase + ".inputs.json").toFile)
      val withInputsAddition = if (inputsFile.isDefined) " and inputs file" else ""
      it should s"successfully validate $versionName workflow: '$validCase'$withInputsAddition" in {
        val wdl = mustExist(validTestCases.resolve(validCase).resolve(validCase + ".wdl").toFile)
        val inputsArgs = inputsFile match {
          case Some(path) => Seq("-i", path.getAbsolutePath)
          case None => Seq.empty[String]
        }

        WomtoolMain.runWomtool(Seq("validate", wdl.getAbsolutePath) ++ inputsArgs) should be(SuccessfulTermination("Success!"))
      }
    }

    listFilesAndFilterDSFile(invalidTestCases) foreach { invalidCase =>
      val inputsFile = ifExists(invalidTestCases.resolve(invalidCase).resolve(invalidCase + ".inputs.json").toFile)
      val withInputsAddition = if (inputsFile.isDefined) " and inputs file" else ""

      it should s"fail to validate $versionName workflow: '$invalidCase'$withInputsAddition" in {
        val wdl = mustExist(invalidTestCases.resolve(invalidCase).resolve(invalidCase + ".wdl").toFile)
        val errorFile = ifExists(invalidTestCases.resolve(invalidCase).resolve("error.txt").toFile).map(f => File(f.getAbsolutePath).contentAsString)
        val inputsArgs = inputsFile match {
          case Some(path) => Seq("-i", path.getAbsolutePath)
          case None => Seq.empty[String]
        }

        WomtoolMain.runWomtool(Seq("validate", wdl.getAbsolutePath) ++ inputsArgs) match {
          case UnsuccessfulTermination(msg) => errorFile match {
            case Some(expectedError) =>
              msg should include(expectedError.trim.replace(s"$${PWD_NAME}", presentWorkingDirectoryName))
            case None => succeed
          }
          case other => fail(s"Expected UnsuccessfulTermination but got $other")
        }
      }
    }
  }


  behavior of "womtool validate with --list-dependencies flag"

  val validationWithImportsTests = File("womtool/src/test/resources/validate-with-imports")
  val validateWithImportsLanguageVersions = Option(validationWithImportsTests.list).toList.flatten
  val userDirectory = sys.props("user.dir")
  val workingDirectory = File(sys.env.getOrElse("CROMWELL_BUILD_ROOT_DIRECTORY", userDirectory)).pathAsString

  it should "test at least one version" in {
    validateWithImportsLanguageVersions.isEmpty should be(false)
  }

  validateWithImportsLanguageVersions.filterNot(f => f.name.contains(".DS")) foreach { versionDirectory =>
    val versionName = versionDirectory.name

    it should s"test at least one $versionName workflow" in {
      versionDirectory.isEmpty should be(false)
    }

    Option(versionDirectory.list).toList.flatten.filterNot(s => s.pathAsString.contains(".DS")) foreach { validCase =>
      val caseName = validCase.name
      val wdlFile = mustExist(versionDirectory.path.resolve(s"../../validate/$versionName/valid/$caseName/$caseName.wdl").toFile)
       wdlFile.getAbsolutePath.split("cromwell/womtool")(0)

      it should s"successfully validate and print the workflow dependencies for $versionName workflow: '$caseName'" in {
        val rawOutput = expectedOutput(versionDirectory, caseName, "expected_imports.txt")
        val importsExpectation = rawOutput.replaceAll("\\{REPLACE_WITH_ROOT_PATH\\}", workingDirectory)

        val res = WomtoolMain.runWomtool(Seq("validate", "-l", wdlFile.getAbsolutePath))
        assert(res.isInstanceOf[SuccessfulTermination])
        val stdout = res.stdout.get
        stdout shouldBe importsExpectation
      }
    }
  }


  private def mustExist(file: java.io.File): java.io.File = if (file.exists) file else fail(s"No such file: ${file.getAbsolutePath}")
  private def ifExists(file: java.io.File): Option[java.io.File] = if (file.exists) Option(file) else None

  // The filterNot(_.contains(".DS")) stuff prevents Mac 'Desktop Services' hidden directories from accidentally being picked up:
  private def listFilesAndFilterDSFile(path: Path): immutable.Seq[String] =  Option(path.toFile.list).toList.flatten.filterNot(_.contains(".DS"))

  private def expectedOutput(versionDirectory: File, caseName: String, outputTextFileName: String): String =
    File(mustExist(versionDirectory.path.resolve(caseName).resolve(outputTextFileName).toFile).getAbsolutePath).contentAsString.trim
}
