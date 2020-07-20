package cromwell.services.womtool

import java.nio.file.Files

import better.files.File
import cromwell.core.{WorkflowOptions, WorkflowSourceFilesCollection, WorkflowSourceFilesWithoutImports}
import cromwell.languages.config.{CromwellLanguages, LanguageConfiguration}
import cromwell.services.womtool.DescriberSpec._
import cromwell.services.womtool.WomtoolServiceMessages.DescribeSuccess
import io.circe.Json
import org.scalatest.{Assertion, FlatSpec, Matchers}

import scala.util.Try

class DescriberSpec extends FlatSpec with Matchers {

  val validationTestCases = File("services/src/test/resources/describe")
  val languageVersions = Option(validationTestCases.list).toList.flatten

  CromwellLanguages.initLanguages(LanguageConfiguration.AllLanguageEntries)

  behavior of "cromwell.services.womtool.Describer"

  it should "test at least one version" in {
    languageVersions.isEmpty should be(false)
  }

  // The filterNot(_.contains(".DS")) stuff prevents Mac 'Desktop Services' hidden directories from accidentally being picked up:
  languageVersions.filterNot(f => f.name.contains(".DS")) foreach { versionDirectory =>
    versionDirectory.path.toFile.listFiles().filterNot(f => f.getName.contains(".DS")) foreach { caseDirectory: java.io.File =>
      it should s" describe ${caseDirectory.getName} (${versionDirectory.name})" in {

        val testCase = DescriberSpec.interpretTestCase(caseDirectory)

        val workflowType = Try(scala.io.Source.fromFile(caseDirectory.toPath.resolve("workflowType").toFile).mkString.stripLineEnd).toOption
        val workflowTypeVersion = Try(scala.io.Source.fromFile(caseDirectory.toPath.resolve("workflowTypeVersion").toFile).mkString.stripLineEnd).toOption

        val interimWsfc = WorkflowSourceFilesWithoutImports(
          workflowSource = None,
          workflowUrl = None,
          workflowRoot = None,
          workflowType = workflowType,
          workflowTypeVersion = workflowTypeVersion,
          inputsJson = "",
          workflowOptions = WorkflowOptions.empty,
          labelsJson = "",
          warnings = Seq.empty
        )

        val wsfc = testCase match {
          case FileAndDescription(file, _) => interimWsfc.copy(workflowSource = Option(file))
          case UrlAndDescription(url, _) => interimWsfc.copy(workflowUrl = Option(url))
        }

        check(wsfc, io.circe.parser.parse(testCase.expectedDescription).getOrElse(???))
      }
    }
  }

  private def check(wsfc: WorkflowSourceFilesCollection, expectedJson: Json): Assertion = {
    import cromwell.services.womtool.models.WorkflowDescription.workflowDescriptionEncoder
    import io.circe.syntax._

    // We test scenarios that would produce DescribeFailure in
    // cromwell.services.womtool.impl.WomtoolServiceInCromwellActorSpec
    // and
    // cromwell.webservice.routes.WomtoolRouteSupportSpec
    // so here we bravely assume that we are able to describe successfully (both valid and invalid workflows)
    //
    // N.B. the `asJson` is highly significant as it exercises the entire serialization module and compares
    // the end product instead of an intermediate case class hierarchy
    Describer.describeWorkflow(wsfc).asInstanceOf[DescribeSuccess].description.asJson shouldBe expectedJson
  }
}

object DescriberSpec {
  sealed trait DescriberSpecTestCase { def expectedDescription: String }
  final case class FileAndDescription(file: String, override val expectedDescription: String) extends DescriberSpecTestCase
  final case class UrlAndDescription(url: String, override val expectedDescription: String) extends DescriberSpecTestCase

  def interpretTestCase(caseDirectory: java.io.File): DescriberSpecTestCase = {
    val description = scala.io.Source.fromFile(caseDirectory.toPath.resolve("description.json").toFile).mkString
    if (Files.exists(caseDirectory.toPath.resolve("workflow.wdl"))) {
      val workflow = scala.io.Source.fromFile(caseDirectory.toPath.resolve("workflow.wdl").toFile).mkString
      FileAndDescription(workflow, description)
    } else if (Files.exists(caseDirectory.toPath.resolve("workflow_url.txt"))) {
      val workflowUrl = scala.io.Source.fromFile(caseDirectory.toPath.resolve("workflow_url.txt").toFile).mkString.trim
      UrlAndDescription(workflowUrl, description)
    } else throw new RuntimeException("Bad test case setup: Expected one of workflow.wdl or workflow_url.txt")


  }
}
