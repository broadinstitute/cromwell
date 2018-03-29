package womtool

import better.files.File
import org.scalatest.{FlatSpec, Matchers}
import womtool.WomtoolMain.{SuccessfulTermination, UnsuccessfulTermination}

class WomtoolValidateDraft2Spec extends FlatSpec with Matchers {
  behavior of "womtool validate for WDL draft 2"

  val wdlDraft2TestCases = File("womtool/src/test/resources/validate/draft2")

  val validTestCases = wdlDraft2TestCases.path.resolve("valid")
  val invalidTestCases = wdlDraft2TestCases.path.resolve("invalid")

  List(validTestCases, invalidTestCases) foreach { path =>
    it should s"be set up for testing '${wdlDraft2TestCases.relativize(path).toString}'" in {
      if (!path.toFile.exists) fail(s"Path doesn't exist: ${path.toAbsolutePath.toString}")
      if (Option(path.toFile.list).toList.flatten.isEmpty) fail(s"No test cases found in: ${path.toAbsolutePath.toString}")
      wdlDraft2TestCases.list.nonEmpty shouldBe true
    }
  }

  Option(validTestCases.toFile.list).toList.flatten foreach { validCase =>
    it should s"successfully validate WDL draft 2: '$validCase'" in {
      val wdl = mustExist(validTestCases.resolve(validCase).resolve(validCase + ".wdl").toFile)
      val inputsArgs = ifExists(invalidTestCases.resolve(validCase).resolve(validCase + ".inputs.json").toFile) match {
        case Some(path) => Seq("-i", path.getAbsolutePath)
        case None => Seq.empty[String]
      }

      WomtoolMain.runWomtool(Seq("validate", wdl.getAbsolutePath) ++ inputsArgs) should be(SuccessfulTermination(""))
    }
  }

  Option(invalidTestCases.toFile.list).toList.flatten foreach { invalidCase =>
    it should s"fail to validate WDL draft 2: '$invalidCase'" in {
      val wdl = mustExist(invalidTestCases.resolve(invalidCase).resolve(invalidCase + ".wdl").toFile)
      val errorFile = ifExists(invalidTestCases.resolve(invalidCase).resolve("error.txt").toFile).map(f => File(f.getAbsolutePath).contentAsString)
      val inputsArgs = ifExists(invalidTestCases.resolve(invalidCase).resolve(invalidCase + ".inputs.json").toFile) match {
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

  private def mustExist(file: java.io.File): java.io.File = if (file.exists) file else fail(s"No such file: ${file.getAbsolutePath}")
  private def ifExists(file: java.io.File): Option[java.io.File] = if (file.exists) Option(file) else None

}
