package womtool

import better.files.File
import org.scalatest.{FlatSpec, Matchers}
import womtool.WomtoolMain.{SuccessfulTermination, UnsuccessfulTermination}

class WomtoolValidateSpec extends FlatSpec with Matchers {

  val validationTestCases = File("womtool/src/test/resources/validate")
  val languageVersions = Option(validationTestCases.list).toList.flatten

  behavior of "womtool validate"


  it should "test at least one version" in {
    languageVersions.isEmpty should be(false)
  }

  languageVersions foreach { versionDirectory =>
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

    Option(validTestCases.toFile.list).toList.flatten foreach { validCase =>
      val inputsFile = ifExists(validTestCases.resolve(validCase).resolve(validCase + ".inputs.json").toFile)
      val withInputsAddition = if (inputsFile.isDefined) " and inputs file" else ""
      it should s"successfully validate $versionName workflow: '$validCase'$withInputsAddition" in {
        val wdl = mustExist(validTestCases.resolve(validCase).resolve(validCase + ".wdl").toFile)
        val inputsArgs = inputsFile match {
          case Some(path) => Seq("-i", path.getAbsolutePath)
          case None => Seq.empty[String]
        }

        WomtoolMain.runWomtool(Seq("validate", wdl.getAbsolutePath) ++ inputsArgs) should be(SuccessfulTermination(""))
      }
    }

    Option(invalidTestCases.toFile.list).toList.flatten foreach { invalidCase =>
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
            case Some(expectedError) => msg should include (expectedError.trim)
            case None => succeed
          }
          case other => fail(s"Expected UnsuccessfulTermination but got $other")
        }
      }
    }
  }

  private def mustExist(file: java.io.File): java.io.File = if (file.exists) file else fail(s"No such file: ${file.getAbsolutePath}")
  private def ifExists(file: java.io.File): Option[java.io.File] = if (file.exists) Option(file) else None

}
