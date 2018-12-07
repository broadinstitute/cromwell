package cromwell.services.womtool.impl

import akka.pattern._
import akka.testkit.TestProbe
import com.typesafe.config.ConfigFactory
import cromwell.core.WorkflowSourceFilesCollection
import cromwell.languages.config.{CromwellLanguages, LanguageConfiguration}
import cromwell.services.ServicesSpec
import cromwell.services.womtool.WomtoolServiceMessages._
import wom.core._

class WomtoolServiceInCromwellActorSpec extends ServicesSpec("Womtool") {

  val womtoolActor = system.actorOf(WomtoolServiceInCromwellActor.props(ConfigFactory.empty(), ConfigFactory.empty(), TestProbe().ref))
  CromwellLanguages.initLanguages(LanguageConfiguration.AllLanguageEntries)

  object TestData {
    val workflowUrlValid = "https://raw.githubusercontent.com/broadinstitute/cromwell/develop/womtool/src/test/resources/validate/wdl_draft3/valid/callable_imports/my_workflow.wdl"
    val workflowUrlNotFound = "https://raw.githubusercontent.com/broadinstitute/cromwell/develop/my_workflow"
    val workflowUrlBadHost = "https://zardoz.zardoz"
    val workflowUrlNotAUrl = "Zardoz"

    val wdlValid =
      s"""
         |task hello {
         |  String addressee
         |  command {
         |    echo "Hello World!"
         |  }
         |  output {
         |    String salutation = read_string(stdout())
         |  }
         |}
         |
         |workflow wf_hello {
         |  call hello
         |}
         |""".stripMargin
    val helloWorldInputs = """{"wf_hello.hello.addressee": "World"}"""
    val bogusInputs = """{"foo.bar": "World"}"""
    val wdlInvalid = "This is not a valid WDL."
  }

  "WomtoolServiceInCromwellActor" should {

    "return valid for a valid workflow" in {

      val wsfc = wsfcConjurer(workflowSource = Some(TestData.wdlValid))

      check(DescribeRequest(wsfc), DescribeSuccess(description = WorkflowDescription(valid = true, errors = List.empty)))
    }

    "return valid for a valid workflow with matching inputs" in {

      val wsfc = wsfcConjurer(workflowSource = Some(TestData.wdlValid), inputsJson = TestData.helloWorldInputs)

      check(DescribeRequest(wsfc), DescribeSuccess(description = WorkflowDescription(valid = true, errors = List.empty)))
    }

    "return invalid for a valid workflow with the wrong inputs" in {

      val wsfc = wsfcConjurer(workflowSource = Some(TestData.wdlValid), inputsJson = TestData.bogusInputs)

      check(DescribeRequest(wsfc), DescribeSuccess(
        description = WorkflowDescription(valid = false, errors = List("Required workflow input 'wf_hello.hello.addressee' not specified"))))
    }

    "return valid for a valid workflow URL" in {

      val wsfc = wsfcConjurer(workflowUrl = Some(TestData.workflowUrlValid))

      check(DescribeRequest(wsfc), DescribeSuccess(description = WorkflowDescription(valid = true, errors = List.empty)))
    }

    "return an error with empty inputs" in {

      val wsfc = wsfcConjurer()

      check(DescribeRequest(wsfc), DescribeFailure("Either workflow source or url has to be supplied"))
    }

    "return an error when both workflow URL and workflow source specified" in {

      val wsfc = wsfcConjurer(workflowSource = Some(TestData.wdlInvalid), workflowUrl = Some(TestData.workflowUrlValid))

      check(DescribeRequest(wsfc), DescribeFailure("Both workflow source and url can't be supplied"))
    }

    "return an error when the workflow URL is a 404" in {

      val wsfc = wsfcConjurer(workflowUrl = Some(TestData.workflowUrlNotFound))

      check(DescribeRequest(wsfc),
        DescribeFailure(
          "Failed to resolve 'https://raw.githubusercontent.com/broadinstitute/cromwell/develop/my_workflow' using resolver: 'http importer (no 'relative-to' origin)' (reason 1 of 1): Failed to download https://raw.githubusercontent.com/broadinstitute/cromwell/develop/my_workflow (reason 1 of 1): 404: Not Found\n"))
    }

    "return an error when the workflow URL's host can't be resolved" in {

      val wsfc = wsfcConjurer(workflowUrl = Some(TestData.workflowUrlBadHost))

      // The error is from the OS network stack and differs between Mac and Linux
      (for {
        result <- (womtoolActor ? DescribeRequest(wsfc)).mapTo[DescribeResult]
        _ = result should (
          be(DescribeFailure("Failed to resolve 'https://zardoz.zardoz' using resolver: 'http importer (no 'relative-to' origin)' (reason 1 of 1): Failed to download https://zardoz.zardoz (reason 1 of 1): HTTP resolver with headers had an unexpected error (zardoz.zardoz: Name or service not known)"))
          or
          be(DescribeFailure("Failed to resolve 'https://zardoz.zardoz' using resolver: 'http importer (no 'relative-to' origin)' (reason 1 of 1): Failed to download https://zardoz.zardoz (reason 1 of 1): HTTP resolver with headers had an unexpected error (zardoz.zardoz: nodename nor servname provided, or not known)"))
        )
      } yield ()).futureValue
    }

    "return an error when the workflow URL is not a URL" in {

      val wsfc = wsfcConjurer(workflowUrl = Some(TestData.workflowUrlNotAUrl))

      // The HTTP resolver has figured out that you have not given it a URL and assumes it's a relative path
      check(DescribeRequest(wsfc), DescribeFailure("Failed to resolve 'Zardoz' using resolver: 'http importer (no 'relative-to' origin)' (reason 1 of 1): Relative path"))
    }

  }

  private def check(request: DescribeRequest, expectedResponse: DescribeResult) = {
    (for {
      result <- (womtoolActor ? request).mapTo[DescribeResult]
      _ = result shouldBe expectedResponse
    } yield ()).futureValue
  }

  private def wsfcConjurer(workflowSource: Option[WorkflowSource] = None,
                           workflowUrl: Option[WorkflowUrl] = None,
                           workflowType: Option[WorkflowType] = None,
                           workflowTypeVersion: Option[WorkflowTypeVersion] = None,
                           inputsJson: WorkflowJson = ""): WorkflowSourceFilesCollection = {
    WorkflowSourceFilesCollection(
      workflowSource = workflowSource,
      workflowUrl = workflowUrl,
      workflowRoot = None,
      workflowType = workflowType,
      workflowTypeVersion = workflowTypeVersion,
      inputsJson = inputsJson,
      workflowOptionsJson = "",
      labelsJson = "",
      importsFile = None,
      workflowOnHold = false,
      warnings = Seq.empty
    )
  }

}
