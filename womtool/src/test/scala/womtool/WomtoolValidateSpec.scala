package womtool

import better.files.File
import cromwell.core.path.DefaultPathBuilder
import org.scalatest.{FlatSpec, Matchers}
import womtool.WomtoolMain.{SuccessfulTermination, UnsuccessfulTermination}

class WomtoolValidateSpec extends FlatSpec with Matchers {

  behavior of "womtool validate"

  private val presentWorkingDirectoryName = DefaultPathBuilder.get(".").toAbsolutePath.name
  val validationTestCases = File("womtool/src/test/resources/validate")
  val languages = Map(
    "wdl" -> Seq(
      "biscayne",
      "wdl_draft2",
      "wdl_draft3"
    ),
    "cwl" -> Seq(
      "cwlV1_0"
    )
  )

  languages foreach { case (language: String, languageVersions: Seq[String]) =>
    languageVersions foreach { languageVersion: String =>
      val versionDirectory = validationTestCases./(language)./(languageVersion)

      testValid(language, languageVersion, versionDirectory)
      testInvalid(language, languageVersion, versionDirectory)
      logIgnored(language, languageVersion, versionDirectory)
    }
  }

  private def logIgnored(language: String, languageVersion: String, versionDirectory: File) = {
    val ignoredTestCases = versionDirectory.path.resolve("ignored")

    Option(ignoredTestCases.toFile.list).toList.filterNot(_.contains(".DS_Store")).flatten foreach { ignoredCase =>
      it should s"run $language (version $languageVersion) test '$ignoredCase'" ignore {}
    }
  }

  private def testInvalid(language: String, languageVersion: String, versionDirectory: File) = {
    val invalidTestCases = versionDirectory.path.resolve("invalid")
    if (invalidTestCases.toFile.list.isEmpty) fail("No invalid test cases in versionDirectory")

    invalidTestCases.toFile.list.filterNot(_.contains(".DS_Store")) foreach { invalidCase: String =>
      val inputsFile = ifExists(invalidTestCases.resolve(invalidCase).resolve(invalidCase + ".inputs.json").toFile)
      val withInputsAddition = if (inputsFile.isDefined) " and inputs file" else ""

      it should s"fail to validate $language (version $languageVersion) workflow: '$invalidCase'$withInputsAddition" in {
        val workflow = mustExist(invalidTestCases.resolve(invalidCase).resolve(s"$invalidCase.$language").toFile)
        val errorFile = ifExists(invalidTestCases.resolve(invalidCase).resolve("error.txt").toFile).map(f => File(f.getAbsolutePath).contentAsString)
        val inputsArgs = inputsFile match {
          case Some(path) => Seq("-i", path.getAbsolutePath)
          case None => Seq.empty[String]
        }
        WomtoolMain.runWomtool(Seq("validate", workflow.getAbsolutePath) ++ inputsArgs) match {
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

  private def testValid(language: String, languageVersion: String, versionDirectory: File) = {
    val validTestCases = versionDirectory.path.resolve("valid")
    if (validTestCases.toFile.list.isEmpty) fail("No valid test cases in versionDirectory")

    validTestCases.toFile.list.filterNot(_.contains(".DS_Store")) foreach { validCase: String =>
      val inputsFile = ifExists(validTestCases.resolve(validCase).resolve(validCase + ".inputs.json").toFile)
      val withInputsAddition = if (inputsFile.isDefined) " and inputs file" else ""

      it should s"successfully validate $language (version $languageVersion) workflow: '$validCase'$withInputsAddition" in {
        val workflow = mustExist(validTestCases.resolve(validCase).resolve(s"$validCase.$language").toFile)
        val inputsArgs = inputsFile match {
          case Some(path) => Seq("-i", path.getAbsolutePath)
          case None => Seq.empty[String]
        }
        WomtoolMain.runWomtool(Seq("validate", workflow.getAbsolutePath) ++ inputsArgs) should be(SuccessfulTermination(""))
      }
    }
  }

  private def mustExist(file: java.io.File): java.io.File = if (file.exists) file else fail(s"No such file: ${file.getAbsolutePath}")
  private def ifExists(file: java.io.File): Option[java.io.File] = if (file.exists) Option(file) else None

}
