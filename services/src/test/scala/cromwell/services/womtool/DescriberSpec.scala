package cromwell.services.womtool

import better.files.File
import cromwell.core.WorkflowSourceFilesCollection
import cromwell.languages.config.{CromwellLanguages, LanguageConfiguration}
import cromwell.services.womtool.WomtoolServiceMessages.DescribeSuccess
import io.circe.Json
import org.scalatest.{Assertion, FlatSpec, Matchers}

class DescriberSpec extends FlatSpec with Matchers {

  val validationTestCases = File("services/src/test/resources/describe")
  val languageVersions = Option(validationTestCases.list).toList.flatten

  CromwellLanguages.initLanguages(LanguageConfiguration.AllLanguageEntries)

  behavior of "womtool validate"

  it should "test at least one version" in {
    languageVersions.isEmpty should be(false)
  }

  // The filterNot(_.contains(".DS")) stuff prevents Mac 'Desktop Services' hidden directories from accidentally being picked up:
  languageVersions.filterNot(f => f.name.contains(".DS")) foreach { versionDirectory =>
    versionDirectory.path.toFile.listFiles().filterNot(f => f.getName.contains(".DS")) foreach { caseDirectory: java.io.File =>
      it should s"successfully describe ${caseDirectory.getName} (${versionDirectory.name})" in {
        val workflow = scala.io.Source.fromFile(caseDirectory.toPath.resolve("workflow.wdl").toFile).mkString
        val description = scala.io.Source.fromFile(caseDirectory.toPath.resolve("description.json").toFile).mkString

        val wsfc = WorkflowSourceFilesCollection(
          workflowSource = Some(workflow),
          workflowUrl = None,
          workflowRoot = None,
          workflowType = None,
          workflowTypeVersion = None,
          inputsJson = "",
          workflowOptionsJson = "",
          labelsJson = "",
          importsFile = None,
          workflowOnHold = false,
          warnings = Seq.empty
        )

        check(wsfc, io.circe.parser.parse(description).getOrElse(???))
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
