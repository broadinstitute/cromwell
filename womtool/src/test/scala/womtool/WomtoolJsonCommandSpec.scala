package womtool

import better.files.File
import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import womtool.WomtoolJsonCommandSpec._
import womtool.WomtoolMain.SuccessfulTermination

object WomtoolJsonCommandSpec {
  final case class TestDefinition(testName: String, commandFormat: Seq[String], expectationFilename: String)

  private def jsonLines(json: String): Set[String] =
    json.linesIterator.toSet[String].map(_.stripSuffix(",")).filterNot(_.forall(_.isWhitespace))
}

// Test for womtool command line commands which output JSON
abstract class WomtoolJsonCommandSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  val womtoolCommand: String
  val testDefinitions: Seq[TestDefinition]

  def defineTests(): Unit = {
    val commandTestCases: File = File(s"womtool/src/test/resources/$womtoolCommand")
    val languageVersions: Seq[File] = Option(commandTestCases.list).toList.flatten

    behavior of s"womtool $womtoolCommand"

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

        // The WDL file is expected to be valid:
        val wdlFile = mustExist(versionDirectory.path.resolve(s"../../validate/$versionName/valid/$caseName/$caseName.wdl").toFile)

        testDefinitions foreach { definition =>

          it should s"validate '${definition.testName}' for $versionName workflow: '$caseName''" in {
            val expectation = expectedJson(versionDirectory, caseName, definition.expectationFilename)
            val fullCommandFormat = definition.commandFormat :+ wdlFile.getAbsolutePath

            WomtoolMain.runWomtool(fullCommandFormat) match {
              case SuccessfulTermination(actualContent) =>
                val actualSet = jsonLines(actualContent)
                val expectedSet = jsonLines(expectation)

                val unexpected = actualSet.diff(expectedSet)
                val ungenerated = expectedSet.diff(actualSet)

                assert(actualSet == expectedSet, s"Received lines: $actualContent${System.lineSeparator}with unexpected values: ${unexpected.mkString("[", ",", "]")}${System.lineSeparator}and missing expected values: ${ungenerated.mkString("[", ",", "]")}")

              case other => fail(s"Expected successful termination but got $other")
            }
          }
        }
      }
    }
  }

  private def expectedJson(versionDirectory: File, caseName: String, jsonName: String): String = {
    File(mustExist(versionDirectory.path.resolve(caseName).resolve(jsonName).toFile).getAbsolutePath).contentAsString
  }

  private def mustExist(file: java.io.File): java.io.File = if (file.exists) file else fail(s"No such file: ${file.getAbsolutePath}")
}
