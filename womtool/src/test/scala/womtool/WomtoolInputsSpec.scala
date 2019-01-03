package womtool

import better.files.File
import org.scalatest.{FlatSpec, Matchers}
import womtool.WomtoolMain.SuccessfulTermination

class WomtoolInputsSpec extends FlatSpec with Matchers {

  val inputsTestCases = File("womtool/src/test/resources/inputs")
  val languageVersions = Option(inputsTestCases.list).toList.flatten

  behavior of "womtool inputs"


  it should "test at least one version" in {
    languageVersions.isEmpty should be(false)
  }

  languageVersions foreach { versionDirectory =>
    val versionName = versionDirectory.name

    it should s"test at least one $versionName workflow" in {
      versionDirectory.isEmpty should be(false)
    }

    Option(versionDirectory.list).toList.flatten foreach { validCase =>
      val caseName = validCase.name

      val wdlFile = mustExist(versionDirectory.path.resolve(s"../../validate/$versionName/valid/$caseName/$caseName.wdl").toFile)

      it should s"make the right 'all inputs' file for $versionName workflow: '$caseName'" in {
        val allInputsExpectation = expectedJson(versionDirectory, caseName, "all.inputs.json")
        WomtoolMain.runWomtool(Seq("inputs", wdlFile.getAbsolutePath)) match {
          case SuccessfulTermination(actualContent) =>
            val actualSet = inputSet(actualContent)
            val expectedSet = inputSet(allInputsExpectation)

            val unexpected = actualSet.diff(expectedSet)
            val ungenerated = expectedSet.diff(actualSet)

            assert(actualSet == expectedSet, s"Received inputs: $actualContent${System.lineSeparator}with unexpected values: ${unexpected.mkString("[", ",", "]")}${System.lineSeparator}and missing expected values: ${ungenerated.mkString("[", ",", "]")}")

          case other => fail(s"Expected successful termination but got $other")
        }
      }

      it should s"make the right 'required inputs' file for $versionName workflow: '$caseName'" in {
        val requiredInputsExpectation = expectedJson(versionDirectory, caseName, "required.inputs.json")
        WomtoolMain.runWomtool(Seq("inputs", "-o", "false", wdlFile.getAbsolutePath)) match {
          case SuccessfulTermination(actualContent) =>
            val actualSet = inputSet(actualContent)
            val expectedSet = inputSet(requiredInputsExpectation)

            val unexpected = actualSet.diff(expectedSet)
            val ungenerated = expectedSet.diff(actualSet)

            assert(inputSet(actualContent) == inputSet(requiredInputsExpectation), s"Received inputs: $actualContent${System.lineSeparator}with unexpected values: ${unexpected.mkString} and missing expected values: ${ungenerated.mkString}")

          case other => fail(s"Expected successful termination but got $other")
        }
      }
    }

  }

  private def inputSet(json: String): Set[String] = json.lines.toSet[String].map(_.stripSuffix(",")).filterNot(_.forall(_.isWhitespace))

  private def expectedJson(versionDirectory: File, caseName: String, jsonName: String): String = {
    File(mustExist(versionDirectory.path.resolve(caseName).resolve(jsonName).toFile).getAbsolutePath).contentAsString
  }

  private def mustExist(file: java.io.File): java.io.File = if (file.exists) file else fail(s"No such file: ${file.getAbsolutePath}")
}
