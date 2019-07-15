package cromwell.services.womtool.impl

import akka.pattern._
import akka.testkit.TestProbe
import com.typesafe.config.ConfigFactory
import cromwell.core.{WorkflowOptions, WorkflowSourceFilesCollection}
import cromwell.languages.config.{CromwellLanguages, LanguageConfiguration}
import cromwell.services.ServicesSpec
import cromwell.services.womtool.WomtoolServiceMessages._
import cromwell.services.womtool.models.{InputDescription, OutputDescription, WorkflowDescription}
import wom.core._
import wom.types.{WomIntegerType, WomStringType}

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
         |version 1.0
         |
         |task hello {
         |  input {
         |    String addressee
         |  }
         |  command <<<
         |    echo "Hello World!"
         |  >>>
         |  output {
         |    String salutation = read_string(stdout())
         |  }
         |}
         |
         |workflow wf_hello {
         |  call hello
         |}
         |""".stripMargin

    val wdlHttpImportValid =
      s"""
         |version 1.0
         |
         |import "https://raw.githubusercontent.com/broadinstitute/cromwell/develop/womtool/src/test/resources/validate/wdl_draft3/valid/callable_imports/my_task.wdl"
         |
         |task hello {
         |  input {
         |    String addressee
         |  }
         |  command <<<
         |    echo "Hello World!"
         |  >>>
         |  output {
         |    String salutation = read_string(stdout())
         |  }
         |}
         |
         |workflow wf_hello {
         |  call hello
         |}
         |""".stripMargin

    val wdlValidNoInputs =
      s"""
         |version 1.0
         |
         |workflow wf_hello {
         |  call hello
         |}
         |
         |task hello {
         |  command <<<
         |    echo "Hello World!"
         |  >>>
         |}
       """.stripMargin

    val wdlValidDraft2NoInputs =
      s"""
         |workflow wf_hello {
         |  call hello
         |}
         |
         |task hello {
         |  command <<<
         |    echo "Hello World!"
         |  >>>
         |}
       """.stripMargin

    val helloWorldInputs = """{"wf_hello.hello.addressee": "World"}"""
    val bogusInputs = """{"foo.bar": "World"}"""
    val emptyInputs = "{}"
    val wdlInvalid = "This is not a valid WDL."

    val successfulDescription = WorkflowDescription(
      valid = true,
      errors = List.empty,
      validWorkflow = true,
      name = "wf_hello",
      inputs = List(InputDescription("wf_hello.hello.addressee", WomStringType, "String", optional = false, default = None)),
      outputs = List.empty,
      images = List.empty,
      submittedDescriptorType = Map(
        "descriptorType" -> "WDL",
        "descriptorTypeVersion" -> "1.0"
      ),
      importedDescriptorTypes = List.empty,
      meta = Map.empty,
      parameterMeta = Map.empty,
      isRunnableWorkflow = true
    )
  }

  "WomtoolServiceInCromwellActor" should {

    "return valid for a valid workflow" in {

      val wsfc = wsfcConjurer(workflowSource = Option(TestData.wdlValid))

      check(
        DescribeRequest(wsfc),
        DescribeSuccess(
          description = TestData.successfulDescription
        )
      )
    }

    "return valid for a valid workflow with HTTP imports" in {
      val wsfc = wsfcConjurer(workflowSource = Option(TestData.wdlHttpImportValid))

      check(
        DescribeRequest(wsfc),
        DescribeSuccess(
          description = TestData.successfulDescription
        )
      )
    }

    "return valid for a valid workflow with matching inputs" in {

      val wsfc = wsfcConjurer(workflowSource = Option(TestData.wdlValid), inputsJson = TestData.helloWorldInputs)

      check(
        DescribeRequest(wsfc),
        DescribeSuccess(
          description = TestData.successfulDescription
        )
      )
    }

    "return valid = false, validWorkflow = true for a valid workflow with the wrong inputs" in {

      val wsfc = wsfcConjurer(workflowSource = Option(TestData.wdlValid), inputsJson = TestData.bogusInputs)

      check(DescribeRequest(wsfc), DescribeSuccess(
        description = WorkflowDescription(valid = false, errors = List("Required workflow input 'wf_hello.hello.addressee' not specified"), validWorkflow = true, name = "wf_hello", inputs = List(InputDescription("wf_hello.hello.addressee", WomStringType, "String", false, None)), submittedDescriptorType = Map("descriptorType" -> "WDL", "descriptorTypeVersion" -> "1.0"), isRunnableWorkflow = true)))
    }

    "return valid = false, validWorkflow = true for a valid inputs-requiring workflow with an empty inputs JSON" in {

      val wsfc = wsfcConjurer(workflowSource = Option(TestData.wdlValid), inputsJson = TestData.emptyInputs)

      check(DescribeRequest(wsfc), DescribeSuccess(
        description = WorkflowDescription(valid = false, errors = List("Required workflow input 'wf_hello.hello.addressee' not specified"), validWorkflow = true, name = "wf_hello", inputs = List(InputDescription("wf_hello.hello.addressee", WomStringType, "String", false, None)), submittedDescriptorType = Map("descriptorType" -> "WDL", "descriptorTypeVersion" -> "1.0"), isRunnableWorkflow = true)))
    }

    "return valid for a valid no-inputs workflow with empty inputs" in {

      val wsfc = wsfcConjurer(workflowSource = Option(TestData.wdlValidNoInputs), inputsJson = TestData.emptyInputs)

      check(
        DescribeRequest(wsfc),
        DescribeSuccess(
          description = WorkflowDescription(
            valid = true,
            errors = List.empty,
            validWorkflow = true,
            name = "wf_hello",
            inputs = List.empty,
            outputs = List.empty,
            images = List.empty,
            submittedDescriptorType = Map(
              "descriptorType" -> "WDL",
              "descriptorTypeVersion" -> "1.0"
            ),
            importedDescriptorTypes = List.empty,
            meta = Map.empty,
            parameterMeta = Map.empty,
            isRunnableWorkflow = true
          )
        )
      )
    }

    "return valid = false, validWorkflow = true for a valid no-inputs workflow with extraneous inputs" in {

      val wsfc = wsfcConjurer(workflowSource = Option(TestData.wdlValidNoInputs), inputsJson = TestData.bogusInputs)

      check(DescribeRequest(wsfc), DescribeSuccess(
        description = WorkflowDescription(valid = false, errors = List("WARNING: Unexpected input provided: foo.bar"), validWorkflow = true, name = "wf_hello", inputs = List.empty, submittedDescriptorType = Map("descriptorType" -> "WDL", "descriptorTypeVersion" -> "1.0"), isRunnableWorkflow = true)))
    }

    // In draft-2 we allow extraneous inputs for legacy reasons - e.g. users put comments in them
    "return valid for a valid no-inputs draft-2 workflow with extraneous inputs" in {

      val wsfc = wsfcConjurer(workflowSource = Option(TestData.wdlValidDraft2NoInputs), inputsJson = TestData.bogusInputs)

      check(
        DescribeRequest(wsfc),
        DescribeSuccess(
          description = WorkflowDescription(
            valid = true,
            errors = List.empty,
            validWorkflow = true,
            name = "wf_hello",
            inputs = List.empty,
            outputs = List.empty,
            images = List.empty,
            submittedDescriptorType = Map(
              "descriptorType" -> "WDL",
              "descriptorTypeVersion" -> "draft-2"
            ),
            importedDescriptorTypes = List.empty,
            meta = Map.empty,
            parameterMeta = Map.empty,
            isRunnableWorkflow = true
          )
        )
      )
    }

    "return valid for a valid workflow URL" in {

      val wsfc = wsfcConjurer(workflowUrl = Option(TestData.workflowUrlValid))

      check(
        DescribeRequest(wsfc),
        DescribeSuccess(
          description = WorkflowDescription(
            valid = true,
            errors = List.empty,
            validWorkflow = true,
            name = "my_workflow",
            inputs = List(InputDescription("i", WomIntegerType, "Int", optional = false, default = None)),
            outputs = List(OutputDescription("o", WomIntegerType, "Int")),
            images = List.empty,
            submittedDescriptorType = Map(
              "descriptorType" -> "WDL",
              "descriptorTypeVersion" -> "1.0"
            ),
            importedDescriptorTypes = List.empty,
            meta = Map.empty,
            parameterMeta = Map.empty,
            isRunnableWorkflow = true
          )
        )
      )
    }

    "return an error for an empty describe request" in {

      val wsfc = wsfcConjurer()

      check(DescribeRequest(wsfc), DescribeFailure("Either workflow source or url has to be supplied"))
    }

    "return an error when both workflow URL and workflow source specified" in {

      val wsfc = wsfcConjurer(workflowSource = Option(TestData.wdlInvalid), workflowUrl = Option(TestData.workflowUrlValid))

      check(DescribeRequest(wsfc), DescribeFailure("Both workflow source and url can't be supplied"))
    }

    "return an error when the workflow URL is a 404" in {

      val wsfc = wsfcConjurer(workflowUrl = Option(TestData.workflowUrlNotFound))

      check(DescribeRequest(wsfc),
        DescribeFailure(
          "Failed to resolve 'https://raw.githubusercontent.com/broadinstitute/cromwell/develop/my_workflow' using resolver: 'http importer (no 'relative-to' origin)' (reason 1 of 1): Failed to download https://raw.githubusercontent.com/broadinstitute/cromwell/develop/my_workflow (reason 1 of 1): 404: Not Found\n"))
    }

    "return an error when the workflow URL's host can't be resolved" in {

      val wsfc = wsfcConjurer(workflowUrl = Option(TestData.workflowUrlBadHost))

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

      val wsfc = wsfcConjurer(workflowUrl = Option(TestData.workflowUrlNotAUrl))

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
      workflowOptions = WorkflowOptions.empty,
      labelsJson = "",
      importsFile = None,
      workflowOnHold = false,
      warnings = Seq.empty
    )
  }

}
