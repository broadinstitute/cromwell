package cromwell.engine.workflow.lifecycle

import akka.actor.Props
import akka.testkit.TestDuration
import com.typesafe.config.ConfigFactory
import cromwell.CromwellTestKitWordSpec
import cromwell.core.labels.{Label, Labels}
import cromwell.core.{WorkflowId, WorkflowOptions, WorkflowSourceFilesWithoutImports}
import cromwell.engine.backend.{BackendConfigurationEntry, CromwellBackends}
import cromwell.engine.workflow.lifecycle.MaterializeWorkflowDescriptorActor.{MaterializeWorkflowDescriptorCommand, MaterializeWorkflowDescriptorFailureResponse, MaterializeWorkflowDescriptorSuccessResponse}
import cromwell.util.SampleWdl.HelloWorld
import org.scalatest.BeforeAndAfter
import org.scalatest.mockito.MockitoSugar
import spray.json.DefaultJsonProtocol._
import spray.json._
import wdl4s.values.WdlString

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
  val wdlSourceWithDocker = HelloWorld.wdlSource(""" runtime { docker: "ubuntu:latest" } """)
  val wdlSourceNoDocker = HelloWorld.wdlSource(""" runtime { } """)
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
        wdlSource = wdlSourceNoDocker,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = validInputsJson,
        workflowOptionsJson = validOptionsFile,
        labelsJson = validCustomLabelsFile)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorSuccessResponse(wfDesc) =>
            wfDesc.id shouldBe workflowId
            wfDesc.name shouldBe "wf_hello"
            wfDesc.namespace.tasks.size shouldBe 1
            wfDesc.knownValues.head shouldBe (("wf_hello.hello.addressee", WdlString("world")))
            wfDesc.backendDescriptor.knownValues.head shouldBe (("wf_hello.hello.addressee", WdlString("world")))
            wfDesc.getWorkflowOption(WorkflowOptions.WriteToCache) shouldBe Option("true")
            wfDesc.getWorkflowOption(WorkflowOptions.ReadFromCache) shouldBe None
            wfDesc.backendDescriptor.customLabels shouldBe Labels("label1" -> "value1", "label2" -> "value2")
            // Default backend assignment is "Local":
            wfDesc.backendAssignments foreach {
              case (call, assignment) if call.task.name.equals("hello") => assignment shouldBe "Local"
              case (call, _) => fail(s"Unexpected call: ${call.task.name}")
            }
            wfDesc.pathBuilders.size shouldBe 1
          case MaterializeWorkflowDescriptorFailureResponse(reason) => fail(s"Materialization failed with $reason")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "assign default runtime attributes" ignore {
      val wdl =
        """
          |task a {
          | command {}
          | runtime { docker: "specified:docker" }
          |}
          |task b { command {} }
          |workflow foo {
          | call a
          | call b
          |}
        """.stripMargin

      val defaultDocker =
        """
          |{
          |  "default_runtime_attributes": {
          |    "docker": "default:docker"
          |  }
          |}
        """.stripMargin
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviorActor, workflowId, importLocalFilesystem = false))
      val sources = WorkflowSourceFilesWithoutImports(
        wdlSource = wdl,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = "{}",
        workflowOptionsJson = defaultDocker,
        labelsJson = validCustomLabelsFile)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorSuccessResponse(wfDesc) =>
            wfDesc.namespace.tasks foreach {
              case task if task.name.equals("a") =>
                task.runtimeAttributes.attrs.size shouldBe 1
                task.runtimeAttributes.attrs.head._2 shouldBe "\"specified:docker\""
              case task if task.name.equals("b") =>
                task.runtimeAttributes.attrs.size shouldBe 1
                task.runtimeAttributes.attrs.head._2 shouldBe "\"default:docker\""
              case task => fail(s"Unexpected task: ${task.name}")
            }
          case MaterializeWorkflowDescriptorFailureResponse(reason) => fail(s"This materialization should not have failed (reason: $reason)")
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
        wdlSource = wdl,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = "{}",
        workflowOptionsJson = "{}",
        labelsJson = validCustomLabelsFile)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, differentDefaultBackendConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorSuccessResponse(wfDesc) =>
            wfDesc.namespace.workflow.taskCalls foreach {
              case call if call.task.name.equals("a") =>
                wfDesc.backendAssignments(call) shouldBe "SpecifiedBackend"
              case call if call.task.name.equals("b") =>
                wfDesc.backendAssignments(call) shouldBe "DefaultBackend"
              case call => fail(s"Unexpected task: ${call.task.name}")
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
        wdlSource = wdl,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = "{}",
        workflowOptionsJson = "{}",
        labelsJson = "{}")
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
        wdlSource = unstructuredFile,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = validInputsJson,
        workflowOptionsJson = validOptionsFile,
        labelsJson = validCustomLabelsFile)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed:\nUnable to load namespace from workflow: ERROR: Finished parsing without consuming all tokens.")
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
        wdlSource = noWorkflowWdl,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = validInputsJson,
        workflowOptionsJson = validOptionsFile,
        labelsJson = validCustomLabelsFile)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed:\nUnable to load namespace from workflow: Namespace does not have a local workflow to run")
          case _: MaterializeWorkflowDescriptorSuccessResponse => fail("This materialization should not have succeeded!")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }


    "reject a taskless WDL source" ignore {
      val noWorkflowWdl =
        """
          |# no task foo { ... } block!!
          |
          |workflow foo {  }
        """.stripMargin
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviorActor, workflowId, importLocalFilesystem = false))
      val badWdlSources = WorkflowSourceFilesWithoutImports(
        wdlSource = noWorkflowWdl,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = validInputsJson,
        workflowOptionsJson = validOptionsFile,
        labelsJson = validCustomLabelsFile)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(badWdlSources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed:\nUnable to load namespace from workflow: Namespace does not have a local workflow to run")
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
        wdlSource = wdlSourceNoDocker,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = validInputsJson,
        workflowOptionsJson = unstructuredFile,
        labelsJson = validCustomLabelsFile)
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
        wdlSource = wdlSourceNoDocker,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = validInputsJson,
        workflowOptionsJson = validOptionsFile,
        labelsJson = unstructuredFile)
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
        wdlSource = wdlSourceNoDocker,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = validInputsJson,
        workflowOptionsJson = validOptionsFile,
        labelsJson = badCustomLabelsFile)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should be(
              s"""Workflow input processing failed:
                |Invalid label: Label1 did not match the regex ${Label.LabelRegexPattern}
                |Invalid label: valu£1 did not match the regex ${Label.LabelRegexPattern}
                |Invalid label: --label2 did not match the regex ${Label.LabelRegexPattern}
                |Invalid label: valuevaluevaluevaluevaluevaluevaluevaluevaluevaluevaluevaluevalue was 65 characters. The maximum is 63""".stripMargin)
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
        wdlSource = wdlSourceNoDocker,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = unstructuredFile,
        workflowOptionsJson = validOptionsFile,
        labelsJson = validCustomLabelsFile)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed:\nWorkflow contains invalid inputs JSON")
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
        wdlSource = wdlSourceNoDocker,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = noInputsJson,
        workflowOptionsJson = validOptionsFile,
        labelsJson = validCustomLabelsFile)
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
        wdlSource = wdl,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = "{}",
        workflowOptionsJson = validOptionsFile,
        labelsJson = validCustomLabelsFile)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed:\nUnable to load namespace from workflow: ERROR: Value for j is not coerceable into a Int")
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
          |
          |  call bar
          |}
        """.stripMargin
      val jsonInput = Map("foo.bad_one" -> "\"gs://this/is/a/bad/gcs/path.txt", "foo.good_one" -> "\"/local/path/is/ok.txt", "foo.bad_two" -> "\"gs://another/bad/gcs/path.txt").toJson.toString
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviorActor, workflowId, importLocalFilesystem = false))
      val sources = WorkflowSourceFilesWithoutImports(
        wdlSource = wdl,
        workflowType = Option("WDL"),
        workflowTypeVersion = None,
        inputsJson = jsonInput,
        workflowOptionsJson = validOptionsFile,
        labelsJson = validCustomLabelsFile)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should equal("Workflow input processing failed:\nInvalid value for File input 'foo.bad_one': \"gs://this/is/a/bad/gcs/path.txt starts with a '\"' " +
              "\nInvalid value for File input 'foo.bad_two': \"gs://another/bad/gcs/path.txt starts with a '\"' ")
          case _: MaterializeWorkflowDescriptorSuccessResponse => fail("This materialization should not have succeeded!")
          case unknown => fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }
  }
}
