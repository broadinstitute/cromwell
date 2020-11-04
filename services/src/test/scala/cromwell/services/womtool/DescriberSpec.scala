package cromwell.services.womtool

import common.assertion.CromwellTimeoutSpec
import cromwell.core.path._
import cromwell.core.{WorkflowOptions, WorkflowSourceFilesCollection, WorkflowSourceFilesWithoutImports}
import cromwell.languages.config.{CromwellLanguages, LanguageConfiguration}
import cromwell.services.womtool.DescriberSpec._
import cromwell.services.womtool.WomtoolServiceMessages.DescribeSuccess
import io.circe.Json
import io.circe.parser._
import org.scalatest.Assertion
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.util.Try

class DescriberSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  private val validationTestCases = DefaultPathBuilder.get("services/src/test/resources/describe")
  private val languageVersions = Option(validationTestCases.list).toList.flatten

  CromwellLanguages.initLanguages(LanguageConfiguration.AllLanguageEntries)

  behavior of "cromwell.services.womtool.Describer"

  it should "test at least one version" in {
    languageVersions.isEmpty should be(false)
  }

  // The filterNot(_.contains(".DS")) stuff prevents Mac 'Desktop Services' hidden directories from accidentally being picked up:
  languageVersions.filterNot(_.name.contains(".DS")) foreach { versionDirectory =>
    versionDirectory.list.toList.filterNot(_.name.contains(".DS")) foreach { caseDirectory =>
      it should s" describe ${caseDirectory.name} (${versionDirectory.name})" in {

        val testCase = DescriberSpec.interpretTestCase(caseDirectory)

        val workflowType = Try(caseDirectory.resolve("workflowType").contentAsString.stripLineEnd).toOption
        val workflowTypeVersion =
          Try(caseDirectory.resolve("workflowTypeVersion").contentAsString.stripLineEnd).toOption

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

        check(wsfc, parse(testCase.expectedDescription).right.get)
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

  def interpretTestCase(caseDirectory: Path): DescriberSpecTestCase = {
    val description = caseDirectory.resolve("description.json").contentAsString
    if (caseDirectory.resolve("workflow.wdl").exists) {
      val workflow = caseDirectory.resolve("workflow.wdl").contentAsString
      FileAndDescription(workflow, description)
    } else if (caseDirectory.resolve("workflow_url.txt").exists) {
      val workflowUrl = caseDirectory.resolve("workflow_url.txt").contentAsString.trim
      UrlAndDescription(workflowUrl, description)
    } else throw new RuntimeException("Bad test case setup: Expected one of workflow.wdl or workflow_url.txt")
  }
}
