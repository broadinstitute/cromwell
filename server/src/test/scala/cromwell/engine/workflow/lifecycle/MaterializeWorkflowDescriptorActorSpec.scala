package cromwell.engine.workflow.lifecycle

import akka.actor.{ActorRef, Props}
import akka.testkit.TestDuration
import cats.data.Validated.{Invalid, Valid}
import com.typesafe.config.ConfigFactory
import cromwell.core.CromwellGraphNode._
import cromwell.core._
import cromwell.core.labels.Labels
import cromwell.engine.backend.{BackendConfigurationEntry, CromwellBackends}
import cromwell.engine.workflow.lifecycle.materialization.MaterializeWorkflowDescriptorActor
import cromwell.engine.workflow.lifecycle.materialization.MaterializeWorkflowDescriptorActor.{
  MaterializeWorkflowDescriptorCommand,
  MaterializeWorkflowDescriptorFailureResponse,
  MaterializeWorkflowDescriptorSuccessResponse
}
import cromwell.languages.util.ImportResolver.GithubImportAuthProvider
import cromwell.services.auth.GithubAuthVendingSupport
import cromwell.util.SampleWdl.HelloWorld
import cromwell.{CromwellTestKitSpec, CromwellTestKitWordSpec}
import org.scalatest.BeforeAndAfter
import spray.json.DefaultJsonProtocol._
import spray.json._
import wom.values.{WomInteger, WomString}

import scala.concurrent.duration._

class MaterializeWorkflowDescriptorActorSpec
    extends CromwellTestKitWordSpec
    with BeforeAndAfter
    with GithubAuthVendingSupport {

  private val ioActor = system.actorOf(SimpleIoActor.props)
  private val workflowId = WorkflowId.randomId()
  private val minimumConf = ConfigFactory
    .parseString(
      """
        |backend {
        |  default = "Local"
        |}
        |""".stripMargin
    )
    .withFallback(CromwellTestKitSpec.DefaultConfig)
  private val differentDefaultBackendConf = ConfigFactory
    .parseString(
      """
        |backend {
        |  default = "DefaultBackend"
        |  // These providers are empty here because the MaterializeWorkflowDescriptorActor won't introspect them:
        |  providers {
        |    DefaultBackend {}
        |    SpecifiedBackend {}
        |  }
        |}
        |""".stripMargin
    )
    .withFallback(CromwellTestKitSpec.DefaultConfig)
  private val azurePrivateWorkflowConfig = ConfigFactory
    .parseString(
      """
        |backend {
        |  default = "Local"
        |}
        |
        |private-workflows {
        |  enabled = true
        |  auth {
        |    azure = true
        |  }
        |}
        |""".stripMargin
    )
    .withFallback(CromwellTestKitSpec.DefaultConfig)
  val unstructuredFile = "fubar badness!"
  private val validOptions = WorkflowOptions.fromJsonString(""" { "write_to_cache": true } """).get
  val validCustomLabelsFile = """{ "label1": "value1", "label2": "value2", "Label1": "valu£1" }"""
  val badCustomLabelsFile =
    """{ "key with characters more than 255-at vero eos et accusamus et iusto odio dignissimos ducimus qui blanditiis praesentium voluptatum deleniti atque corrupti quos dolores et quas molestias excepturi sint occaecati cupiditate non provident, similique sunt in culpas": "value with characters more than 255-at vero eos et accusamus et iusto odio dignissimos ducimus qui blanditiis praesentium voluptatum deleniti atque corrupti quos dolores et quas molestias excepturi sint occaecati cupiditate non provident, similique sunt in culpa" }"""

  private val validInputsJson = HelloWorld.rawInputs.toJson.toString()
  private val workflowSourceNoDocker = HelloWorld.workflowSource(""" runtime { } """)
  private val Timeout = 10.second.dilated
  private val NoBehaviorActor = system.actorOf(Props.empty)
  val callCachingEnabled = true
  val invalidateBadCacheResults = true

  private val validMemoryRetryOptions1 =
    WorkflowOptions.fromJsonString(""" { "memory_retry_multiplier": 1.0 } """).get
  private val validMemoryRetryOptions2 =
    WorkflowOptions.fromJsonString(""" { "memory_retry_multiplier": 99.0 } """).get
  private val validMemoryRetryOptions3 =
    WorkflowOptions.fromJsonString(""" { "memory_retry_multiplier": 12.34 } """).get
  private val invalidMemoryRetryOptions1 =
    WorkflowOptions.fromJsonString(""" { "memory_retry_multiplier": 0.9 } """).get
  private val invalidMemoryRetryOptions2 =
    WorkflowOptions.fromJsonString(""" { "memory_retry_multiplier": 99.1 } """).get
  private val invalidMemoryRetryOptions3 =
    WorkflowOptions.fromJsonString(""" { "memory_retry_multiplier": -1.1 } """).get
  private val invalidMemoryRetryOptions4 =
    WorkflowOptions.fromJsonString(""" { "memory_retry_multiplier": "invalid value" } """).get
  private val invalidMemoryRetryOptions5 =
    WorkflowOptions.fromJsonString(""" { "memory_retry_multiplier": true } """).get

  before {}

  after {
    system.stop(NoBehaviorActor)
    system.stop(ioActor)
  }

  private val fooHogGroup = HogGroup("foo")

  override def serviceRegistryActor: ActorRef = NoBehaviorActor

  "MaterializeWorkflowDescriptorActor" should {
    "accept valid WDL, inputs and options files" in {
      val materializeWfActor = system.actorOf(
        MaterializeWorkflowDescriptorActor.props(NoBehaviorActor,
                                                 workflowId,
                                                 importLocalFilesystem = false,
                                                 ioActorProxy = ioActor,
                                                 hogGroup = fooHogGroup
        )
      )
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = Option(workflowSourceNoDocker),
        workflowUrl = None,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = validInputsJson,
        workflowOptions = validOptions,
        labelsJson = validCustomLabelsFile,
        warnings = Vector.empty,
        requestedWorkflowId = None
      )
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources,
                                                                minimumConf,
                                                                callCachingEnabled,
                                                                invalidateBadCacheResults
      )

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorSuccessResponse(wfDesc) =>
            wfDesc.id shouldBe workflowId
            wfDesc.name shouldBe "wf_hello"
            wfDesc.callable.taskCallNodes.size shouldBe 1
            wfDesc.knownValues.head._1.fullyQualifiedName shouldBe "wf_hello.hello.addressee"
            wfDesc.knownValues.head._2 shouldBe WomString("world")
            wfDesc.getWorkflowOption(WorkflowOptions.WriteToCache) shouldBe Option("true")
            wfDesc.getWorkflowOption(WorkflowOptions.ReadFromCache) shouldBe None
            wfDesc.backendDescriptor.customLabels shouldBe
              Labels("Label1" -> "valu£1", "label1" -> "value1", "label2" -> "value2")
            // Default backend assignment is "Local":
            wfDesc.backendAssignments foreach {
              case (call, assignment) if call.callable.name.equals("hello") => assignment shouldBe "Local"
              case (call, _) => fail(s"Unexpected call: ${call.callable.name}")
            }
            wfDesc.pathBuilders.size shouldBe 2 // one each for the local and http filesystems
          case MaterializeWorkflowDescriptorFailureResponse(reason) => fail(s"Materialization failed with $reason")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "accept valid workflowUrl" in {
      val workflowUrl = Option(
        "https://raw.githubusercontent.com/broadinstitute/cromwell/develop/womtool/src/test/resources/validate/wdl_draft3/valid/callable_imports/my_workflow.wdl"
      )
      val inputs = Map("my_workflow.i" -> 5)
      val materializeWfActor = system.actorOf(
        MaterializeWorkflowDescriptorActor.props(NoBehaviorActor,
                                                 workflowId,
                                                 importLocalFilesystem = false,
                                                 ioActorProxy = ioActor,
                                                 hogGroup = fooHogGroup
        )
      )
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = None,
        workflowUrl = workflowUrl,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = inputs.toJson.toString(),
        workflowOptions = WorkflowOptions.empty,
        labelsJson = "{}",
        warnings = Vector.empty,
        requestedWorkflowId = None
      )
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources,
                                                                minimumConf,
                                                                callCachingEnabled,
                                                                invalidateBadCacheResults
      )

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorSuccessResponse(wfDesc) =>
            wfDesc.id shouldBe workflowId
            wfDesc.name shouldBe "my_workflow"
            wfDesc.callable.taskCallNodes.size shouldBe 0
            wfDesc.knownValues.head._1.fullyQualifiedName shouldBe "i"
            wfDesc.knownValues.head._2 shouldBe WomInteger(5)
            wfDesc.getWorkflowOption(WorkflowOptions.WriteToCache) shouldBe None
            wfDesc.getWorkflowOption(WorkflowOptions.ReadFromCache) shouldBe None
            wfDesc.pathBuilders.size shouldBe 2 // // one each for the local and http filesystems
          case MaterializeWorkflowDescriptorFailureResponse(reason) => fail(s"Materialization failed with $reason")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "assign backends based on runtime attributes" in {
      val wdl =
        """
          |task a {
          | command {}
          | runtime { backend: "SpecifiedBackend" }
          |}
          |task b { command {} }
          |workflow foo {
          | call a
          | call b
          |}
        """.stripMargin

      // The backends need to be in CromwellBackends in order for MaterializeWorkflowDescriptorActor to accept them:
      val placeholdingActorFactoryClass = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      val fauxBackendEntries = List(
        BackendConfigurationEntry("SpecifiedBackend", placeholdingActorFactoryClass, ConfigFactory.parseString("")),
        BackendConfigurationEntry("DefaultBackend", placeholdingActorFactoryClass, ConfigFactory.parseString(""))
      )
      val cromwellBackends = CromwellBackends(fauxBackendEntries)

      // Run the test:
      val materializeWfActor = system.actorOf(
        MaterializeWorkflowDescriptorActor.props(NoBehaviorActor,
                                                 workflowId,
                                                 cromwellBackends,
                                                 importLocalFilesystem = false,
                                                 ioActorProxy = ioActor,
                                                 hogGroup = fooHogGroup
        )
      )
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = Option(wdl),
        workflowUrl = None,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = "{}",
        workflowOptions = WorkflowOptions.empty,
        labelsJson = validCustomLabelsFile,
        warnings = Vector.empty,
        requestedWorkflowId = None
      )
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources,
                                                                differentDefaultBackendConf,
                                                                callCachingEnabled,
                                                                invalidateBadCacheResults
      )

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorSuccessResponse(wfDesc) =>
            wfDesc.callable.taskCallNodes foreach {
              case call if call.callable.name.equals("a") =>
                wfDesc.backendAssignments(call) shouldBe "SpecifiedBackend"
              case call if call.callable.name.equals("b") =>
                wfDesc.backendAssignments(call) shouldBe "DefaultBackend"
              case call => fail(s"Unexpected task: ${call.callable.name}")
            }
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            fail(s"Materialization unexpectedly failed ($reason)")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "reject backend assignment to non-existent backends" in {
      val wdl =
        """
          |task a {
          | command {}
          | runtime { backend: "NoSuchBackend" }
          |}
          |workflow foo {
          | call a
          |}
        """.stripMargin

      val materializeWfActor = system.actorOf(
        MaterializeWorkflowDescriptorActor.props(NoBehaviorActor,
                                                 workflowId,
                                                 importLocalFilesystem = false,
                                                 ioActorProxy = ioActor,
                                                 hogGroup = fooHogGroup
        )
      )
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = Option(wdl),
        workflowUrl = None,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = "{}",
        workflowOptions = WorkflowOptions.empty,
        labelsJson = "{}",
        warnings = Vector.empty,
        requestedWorkflowId = None
      )
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources,
                                                                differentDefaultBackendConf,
                                                                callCachingEnabled,
                                                                invalidateBadCacheResults
      )

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            if (
              !reason.getMessage.contains(
                "Backend for call foo.a ('NoSuchBackend') not registered in configuration file"
              )
            )
              fail(s"Unexpected failure message from MaterializeWorkflowDescriptorActor: ${reason.getMessage}")
          case _: MaterializeWorkflowDescriptorSuccessResponse =>
            fail("This materialization should not have succeeded!")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "reject an invalid WDL source" in {
      val materializeWfActor = system.actorOf(
        MaterializeWorkflowDescriptorActor.props(NoBehaviorActor,
                                                 workflowId,
                                                 importLocalFilesystem = false,
                                                 ioActorProxy = ioActor,
                                                 hogGroup = fooHogGroup
        )
      )
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = Option(unstructuredFile),
        workflowUrl = None,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = validInputsJson,
        workflowOptions = validOptions,
        labelsJson = validCustomLabelsFile,
        warnings = Vector.empty,
        requestedWorkflowId = None
      )
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources,
                                                                minimumConf,
                                                                callCachingEnabled,
                                                                invalidateBadCacheResults
      )

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith(
              "Workflow input processing failed:\nERROR: Finished parsing without consuming all tokens."
            )
          case _: MaterializeWorkflowDescriptorSuccessResponse =>
            fail("This materialization should not have succeeded!")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "reject workflow with invalid URL" in {
      val workflowUrl = Option("https://raw.githubusercontent.com/broadinstitute/cromwell/develop/my_workflow")
      val inputs = Map("my_workflow.i" -> 5)
      val materializeWfActor = system.actorOf(
        MaterializeWorkflowDescriptorActor.props(NoBehaviorActor,
                                                 workflowId,
                                                 importLocalFilesystem = false,
                                                 ioActorProxy = ioActor,
                                                 hogGroup = fooHogGroup
        )
      )
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = None,
        workflowUrl = workflowUrl,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = inputs.toJson.toString(),
        workflowOptions = WorkflowOptions.empty,
        labelsJson = "{}",
        warnings = Vector.empty,
        requestedWorkflowId = None
      )
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources,
                                                                minimumConf,
                                                                callCachingEnabled,
                                                                invalidateBadCacheResults
      )

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith(
              """Workflow input processing failed:
                |Failed to resolve 'https://raw.githubusercontent.com/broadinstitute/cromwell/develop/my_workflow' using resolver: 'http importer (no 'relative-to' origin)' (reason 1 of 1): Failed to download https://raw.githubusercontent.com/broadinstitute/cromwell/develop/my_workflow (reason 1 of 1): 404: Not Found""".stripMargin
            )
          case _: MaterializeWorkflowDescriptorSuccessResponse =>
            fail("This materialization should not have succeeded!")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "reject a workflowless WDL source" in {
      val noWorkflowWdl =
        """
          |task hello { command <<< echo blah >>> }
          |
          |# no workflow foo { ... } block!!
        """.stripMargin
      val materializeWfActor = system.actorOf(
        MaterializeWorkflowDescriptorActor.props(NoBehaviorActor,
                                                 workflowId,
                                                 importLocalFilesystem = false,
                                                 ioActorProxy = ioActor,
                                                 hogGroup = fooHogGroup
        )
      )
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = Option(noWorkflowWdl),
        workflowUrl = None,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = validInputsJson,
        workflowOptions = validOptions,
        labelsJson = validCustomLabelsFile,
        warnings = Vector.empty,
        requestedWorkflowId = None
      )
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources,
                                                                minimumConf,
                                                                callCachingEnabled,
                                                                invalidateBadCacheResults
      )

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith(
              "Workflow input processing failed:\nNamespace does not have a local workflow to run"
            )
          case _: MaterializeWorkflowDescriptorSuccessResponse =>
            fail("This materialization should not have succeeded!")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "reject an invalid workflow inputs file" in {
      val materializeWfActor = system.actorOf(
        MaterializeWorkflowDescriptorActor.props(NoBehaviorActor,
                                                 workflowId,
                                                 importLocalFilesystem = false,
                                                 ioActorProxy = ioActor,
                                                 hogGroup = fooHogGroup
        )
      )
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = Option(workflowSourceNoDocker),
        workflowUrl = None,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = unstructuredFile,
        workflowOptions = validOptions,
        labelsJson = validCustomLabelsFile,
        warnings = Vector.empty,
        requestedWorkflowId = None
      )
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources,
                                                                minimumConf,
                                                                callCachingEnabled,
                                                                invalidateBadCacheResults
      )

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed:\n")
          case _: MaterializeWorkflowDescriptorSuccessResponse =>
            fail("This materialization should not have succeeded!")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "reject requests if any required inputs are missing" in {
      val materializeWfActor = system.actorOf(
        MaterializeWorkflowDescriptorActor.props(NoBehaviorActor,
                                                 workflowId,
                                                 importLocalFilesystem = false,
                                                 ioActorProxy = ioActor,
                                                 hogGroup = fooHogGroup
        )
      )
      val noInputsJson = "{}"
      val badOptionsSources = WorkflowSourceFilesWithoutImports(
        workflowSource = Option(workflowSourceNoDocker),
        workflowUrl = None,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = noInputsJson,
        workflowOptions = validOptions,
        labelsJson = validCustomLabelsFile,
        warnings = Vector.empty,
        requestedWorkflowId = None
      )
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(badOptionsSources,
                                                                minimumConf,
                                                                callCachingEnabled,
                                                                invalidateBadCacheResults
      )

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith(
              "Workflow input processing failed:\nRequired workflow input 'wf_hello.hello.addressee' not specified"
            )
          case _: MaterializeWorkflowDescriptorSuccessResponse =>
            fail("This materialization should not have succeeded!")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "handle coercion failures gracefully" in {
      val wdl =
        """
          |task bar { command { echo foobar } }
          |workflow foo {
          |  Int j = "twenty-seven point five recurring"
          |  call bar
          |}
        """.stripMargin
      val materializeWfActor = system.actorOf(
        MaterializeWorkflowDescriptorActor.props(NoBehaviorActor,
                                                 workflowId,
                                                 importLocalFilesystem = false,
                                                 ioActorProxy = ioActor,
                                                 hogGroup = fooHogGroup
        )
      )
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = Option(wdl),
        workflowUrl = None,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = "{}",
        workflowOptions = validOptions,
        labelsJson = validCustomLabelsFile,
        warnings = Vector.empty,
        requestedWorkflowId = None
      )
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources,
                                                                minimumConf,
                                                                callCachingEnabled,
                                                                invalidateBadCacheResults
      )

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith(
              "Workflow input processing failed:\nERROR: Value 'j' is declared as a 'Int' but the expression evaluates to 'String'"
            )
          case _: MaterializeWorkflowDescriptorSuccessResponse =>
            fail("This materialization should not have succeeded!")
          case unknown => fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "identify all malformed input file names in an input json" in {
      val wdl =
        """
          |task bar { command { echo foobar } }
          |workflow foo {
          |  File bad_one
          |  File good_one
          |  File bad_two
          |  File bad_three
          |
          |  call bar
          |}
        """.stripMargin
      val jsonInput = Map(
        "foo.bad_one" -> "\"gs://this/is/a/bad/gcs/path.txt",
        "foo.good_one" -> "\"/local/path/is/ok.txt",
        "foo.bad_two" -> "\"gs://another/bad/gcs/path.txt",
        "foo.bad_three" -> ""
      ).toJson.toString
      val materializeWfActor = system.actorOf(
        MaterializeWorkflowDescriptorActor.props(NoBehaviorActor,
                                                 workflowId,
                                                 importLocalFilesystem = false,
                                                 ioActorProxy = ioActor,
                                                 hogGroup = fooHogGroup
        )
      )
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = Option(wdl),
        workflowUrl = None,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = jsonInput,
        workflowOptions = validOptions,
        labelsJson = validCustomLabelsFile,
        warnings = Vector.empty,
        requestedWorkflowId = None
      )
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources,
                                                                minimumConf,
                                                                callCachingEnabled,
                                                                invalidateBadCacheResults
      )

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed:\n")
            reason.getMessage should include(
              "Invalid value for File input 'foo.bad_one': \"gs://this/is/a/bad/gcs/path.txt starts with a '\"'"
            )
            reason.getMessage should include(
              "Invalid value for File input 'foo.bad_two': \"gs://another/bad/gcs/path.txt starts with a '\"'"
            )
            reason.getMessage should include("Invalid value for File input 'foo.bad_three': empty value")
          case _: MaterializeWorkflowDescriptorSuccessResponse =>
            fail("This materialization should not have succeeded!")
          case unknown => fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "accept valid memory_retry_multiplier" in {
      List(validMemoryRetryOptions1, validMemoryRetryOptions2, validMemoryRetryOptions3, validOptions) map { options =>
        MaterializeWorkflowDescriptorActor.validateMemoryRetryMultiplier(options) match {
          case Valid(_) => // good!
          case Invalid(_) => fail(s"memory_retry_multiplier validation for $options failed but should have passed!")
        }
      }
    }

    "reject invalid memory_retry_multiplier" in {
      List(invalidMemoryRetryOptions1, invalidMemoryRetryOptions2, invalidMemoryRetryOptions3) map { options =>
        MaterializeWorkflowDescriptorActor.validateMemoryRetryMultiplier(options) match {
          case Invalid(errorsList) =>
            errorsList.head should be(
              "Workflow option 'memory_retry_multiplier' is invalid. " +
                "It should be in the range 1.0 ≤ n ≤ 99.0"
            )
          case Valid(_) => fail(s"memory_retry_multiplier validation for $options succeeded but should have failed!")
        }
      }

      List(invalidMemoryRetryOptions4, invalidMemoryRetryOptions5) map { options =>
        MaterializeWorkflowDescriptorActor.validateMemoryRetryMultiplier(options) match {
          case Invalid(errorsList) =>
            errorsList.head should startWith(
              s"Workflow option 'memory_retry_multiplier' is invalid. " +
                "It should be of type Double and in the range 1.0 ≤ n ≤ 99.0. Error: NumberFormatException:"
            )
          case Valid(_) => fail(s"memory_retry_multiplier validation for $options succeeded but should have failed!")
        }
      }
    }

    "fail materialization if memory_retry_multiplier is invalid" in {
      val materializeWfActor = system.actorOf(
        MaterializeWorkflowDescriptorActor.props(NoBehaviorActor,
                                                 workflowId,
                                                 importLocalFilesystem = false,
                                                 ioActorProxy = ioActor,
                                                 hogGroup = fooHogGroup
        )
      )
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = Option(workflowSourceNoDocker),
        workflowUrl = None,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = validInputsJson,
        workflowOptions = invalidMemoryRetryOptions1,
        labelsJson = validCustomLabelsFile,
        warnings = Vector.empty,
        requestedWorkflowId = None
      )
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources,
                                                                minimumConf,
                                                                callCachingEnabled,
                                                                invalidateBadCacheResults
      )

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed:\n")
            reason.getMessage should include(
              "Workflow option 'memory_retry_multiplier' is invalid. It should be in the range 1.0 ≤ n ≤ 99.0"
            )
          case _: MaterializeWorkflowDescriptorSuccessResponse =>
            fail("This materialization should not have succeeded!")
          case unknown => fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "return Github import auth provider when private workflows in Azure is enabled" in {
      MaterializeWorkflowDescriptorActor.getImportAuthProviders(azurePrivateWorkflowConfig, importAuthProvider) match {
        case Valid(providers) =>
          providers.size shouldBe 1
          providers.head.isInstanceOf[GithubImportAuthProvider] shouldBe true
          providers.head.validHosts shouldBe List("github.com", "githubusercontent.com", "raw.githubusercontent.com")
        case Invalid(e) => fail(s"Unexpected failure: $e")
      }
    }

    "return no import auth provider when private workflows is disabled" in {
      MaterializeWorkflowDescriptorActor.getImportAuthProviders(minimumConf, importAuthProvider) match {
        case Valid(providers) => providers.size shouldBe 0
        case Invalid(e) => fail(s"Unexpected failure: $e")
      }
    }
  }
}
