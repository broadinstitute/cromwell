package cromwell.engine.workflow.lifecycle

import akka.actor.Props
import akka.testkit.TestDuration
import com.typesafe.config.ConfigFactory
import cromwell.CromwellTestKitWordSpec
import cromwell.core.CromwellGraphNode._
import cromwell.core.labels.{Label, Labels}
import cromwell.core.{WorkflowId, WorkflowOptions, WorkflowSourceFilesWithoutImports}
import cromwell.engine.backend.{BackendConfigurationEntry, CromwellBackends}
import cromwell.engine.workflow.lifecycle.materialization.MaterializeWorkflowDescriptorActor
import cromwell.engine.workflow.lifecycle.materialization.MaterializeWorkflowDescriptorActor.{MaterializeWorkflowDescriptorCommand, MaterializeWorkflowDescriptorFailureResponse, MaterializeWorkflowDescriptorSuccessResponse}
import cromwell.util.SampleWdl.HelloWorld
import org.scalatest.BeforeAndAfter
import org.scalatest.mockito.MockitoSugar
import spray.json.DefaultJsonProtocol._
import spray.json._
import wom.values.WomString

import scala.concurrent.duration._

class MaterializeWorkflowDescriptorActorSpec extends CromwellTestKitWordSpec with BeforeAndAfter with MockitoSugar {

  val workflowId = WorkflowId.randomId()
  val minimumConf = ConfigFactory.parseString(
    """
      |backend {
      |  default = "Local"
      |}
    """.stripMargin)
  val differentDefaultBackendConf = ConfigFactory.parseString(
    """
      |backend {
      |  default = "DefaultBackend"
      |  // These providers are empty here because the MaterializeWorkflowDescriptorActor won't introspect them:
      |  providers {
      |    DefaultBackend {}
      |    SpecifiedBackend {}
      |  }
      |}
    """.stripMargin)
  val unstructuredFile = "fubar badness!"
  val validOptionsFile =""" { "write_to_cache": "true" } """
  val validCustomLabelsFile="""{ "label1": "value1", "label2": "value2" }"""
  val badCustomLabelsFile="""{ "Label1": "valu£1", "--label2": "valuevaluevaluevaluevaluevaluevaluevaluevaluevaluevaluevaluevalue" }"""

  val validInputsJson = HelloWorld.rawInputs.toJson.toString()
  val workflowSourceWithDocker = HelloWorld.workflowSource(""" runtime { docker: "ubuntu:latest" } """)
  val workflowSourceNoDocker = HelloWorld.workflowSource(""" runtime { } """)
  val Timeout = 10.second.dilated
  val NoBehaviorActor = system.actorOf(Props.empty)

  before {
  }

  after {
    system.stop(NoBehaviorActor)
  }

  "MaterializeWorkflowDescriptorActor" should {
    "accept valid WDL, inputs and options files" in {
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviorActor, workflowId, importLocalFilesystem = false))
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = workflowSourceNoDocker,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = validInputsJson,
        workflowOptionsJson = validOptionsFile,
        labelsJson = validCustomLabelsFile,
        warnings = Vector.empty)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

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
            wfDesc.backendDescriptor.customLabels shouldBe Labels("label1" -> "value1", "label2" -> "value2")
            // Default backend assignment is "Local":
            wfDesc.backendAssignments foreach {
              case (call, assignment) if call.callable.name.equals("hello") => assignment shouldBe "Local"
              case (call, _) => fail(s"Unexpected call: ${call.callable.name}")
            }
            wfDesc.pathBuilders.size shouldBe 1
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
        BackendConfigurationEntry("DefaultBackend", placeholdingActorFactoryClass, ConfigFactory.parseString("")))
      val cromwellBackends = CromwellBackends(fauxBackendEntries)

      // Run the test:
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviorActor, workflowId, cromwellBackends, importLocalFilesystem = false))
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = wdl,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = "{}",
        workflowOptionsJson = "{}",
        labelsJson = validCustomLabelsFile,
        warnings = Vector.empty)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, differentDefaultBackendConf)

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
          case MaterializeWorkflowDescriptorFailureResponse(reason) => fail(s"Materialization unexpectedly failed ($reason)")
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

      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviorActor, workflowId, importLocalFilesystem = false))
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = wdl,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = "{}",
        workflowOptionsJson = "{}",
        labelsJson = "{}",
        warnings = Vector.empty)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, differentDefaultBackendConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            if (!reason.getMessage.contains("Backend for call foo.a ('NoSuchBackend') not registered in configuration file"))
              fail(s"Unexpected failure message from MaterializeWorkflowDescriptorActor: ${reason.getMessage}")
          case _: MaterializeWorkflowDescriptorSuccessResponse => fail("This materialization should not have succeeded!")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "reject an invalid WDL source" in {
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviorActor, workflowId, importLocalFilesystem = false))
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = unstructuredFile,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = validInputsJson,
        workflowOptionsJson = validOptionsFile,
        labelsJson = validCustomLabelsFile,
        warnings = Vector.empty)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed:\nERROR: Finished parsing without consuming all tokens.")
          case _: MaterializeWorkflowDescriptorSuccessResponse => fail("This materialization should not have succeeded!")
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
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviorActor, workflowId, importLocalFilesystem = false))
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = noWorkflowWdl,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = validInputsJson,
        workflowOptionsJson = validOptionsFile,
        labelsJson = validCustomLabelsFile,
        warnings = Vector.empty)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed:\nNamespace does not have a local workflow to run")
          case _: MaterializeWorkflowDescriptorSuccessResponse => fail("This materialization should not have succeeded!")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "reject an invalid options file" in {
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviorActor, workflowId, importLocalFilesystem = false))
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = workflowSourceNoDocker,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = validInputsJson,
        workflowOptionsJson = unstructuredFile,
        labelsJson = validCustomLabelsFile,
        warnings = Vector.empty)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed:\nWorkflow contains invalid options JSON")
          case _: MaterializeWorkflowDescriptorSuccessResponse => fail("This materialization should not have succeeded!")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "reject an unstructured labels file" in {
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviorActor, workflowId, importLocalFilesystem = false))
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = workflowSourceNoDocker,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = validInputsJson,
        workflowOptionsJson = validOptionsFile,
        labelsJson = unstructuredFile,
        warnings = Vector.empty)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith(
              """Workflow input processing failed:
                |Workflow contains invalid labels JSON: Unexpected character 'u'""".stripMargin)
          case _: MaterializeWorkflowDescriptorSuccessResponse => fail("This materialization should not have succeeded!")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "reject invalid labels" in {
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviorActor, workflowId, importLocalFilesystem = false))
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = workflowSourceNoDocker,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = validInputsJson,
        workflowOptionsJson = validOptionsFile,
        labelsJson = badCustomLabelsFile,
        warnings = Vector.empty)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            val expectedMessage =
              s"""Workflow input processing failed:
                  |Invalid label: `Label1` did not match the regex ${Label.LabelKeyRegex}.
                  |Invalid label: `valu£1` did not match the regex ${Label.LabelValueRegex}.
                  |Invalid label: `--label2` did not match the regex ${Label.LabelKeyRegex}.
                  |Invalid label: `valuevaluevaluevaluevaluevaluevaluevaluevaluevaluevaluevaluevalue` is 65 characters. The maximum is 63.""".stripMargin
            reason.getMessage shouldBe expectedMessage
          case _: MaterializeWorkflowDescriptorSuccessResponse => fail("This materialization should not have succeeded!")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "reject an invalid workflow inputs file" in {
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviorActor, workflowId, importLocalFilesystem = false))
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = workflowSourceNoDocker,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = unstructuredFile,
        workflowOptionsJson = validOptionsFile,
        labelsJson = validCustomLabelsFile,
        warnings = Vector.empty)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed:\n")
          case _: MaterializeWorkflowDescriptorSuccessResponse => fail("This materialization should not have succeeded!")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "reject requests if any required inputs are missing" in {
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviorActor, workflowId, importLocalFilesystem = false))
      val noInputsJson = "{}"
      val badOptionsSources = WorkflowSourceFilesWithoutImports(
        workflowSource = workflowSourceNoDocker,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = noInputsJson,
        workflowOptionsJson = validOptionsFile,
        labelsJson = validCustomLabelsFile,
        warnings = Vector.empty)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(badOptionsSources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed:\nRequired workflow input 'wf_hello.hello.addressee' not specified")
          case _: MaterializeWorkflowDescriptorSuccessResponse => fail("This materialization should not have succeeded!")
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
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviorActor, workflowId, importLocalFilesystem = false))
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = wdl,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = "{}",
        workflowOptionsJson = validOptionsFile,
        labelsJson = validCustomLabelsFile,
        warnings = Vector.empty)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed:\nERROR: Value 'j' is declared as a 'Int' but the expression evaluates to 'String'")
          case _: MaterializeWorkflowDescriptorSuccessResponse => fail("This materialization should not have succeeded!")
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
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviorActor, workflowId, importLocalFilesystem = false))
      val sources = WorkflowSourceFilesWithoutImports(
        workflowSource = wdl,
        workflowRoot = None,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = jsonInput,
        workflowOptionsJson = validOptionsFile,
        labelsJson = validCustomLabelsFile,
        warnings = Vector.empty)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed:\n")
            reason.getMessage should include("Invalid value for File input 'foo.bad_one': \"gs://this/is/a/bad/gcs/path.txt starts with a '\"'")
            reason.getMessage should include("Invalid value for File input 'foo.bad_two': \"gs://another/bad/gcs/path.txt starts with a '\"'")
            reason.getMessage should include("Invalid value for File input 'foo.bad_three': empty value")
          case _: MaterializeWorkflowDescriptorSuccessResponse => fail("This materialization should not have succeeded!")
          case unknown => fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }
  }
}
