package cromwell.engine.workflow

import akka.actor.Props
import akka.testkit.TestDuration
import com.typesafe.config.ConfigFactory
import cromwell.core.Tags._
import cromwell.CromwellTestkitSpec
import cromwell.core.{WorkflowId, WorkflowOptions, WorkflowSourceFiles}
import cromwell.engine.backend.{BackendConfigurationEntry, CromwellBackends}
import cromwell.engine.workflow.lifecycle.MaterializeWorkflowDescriptorActor
import cromwell.engine.workflow.lifecycle.MaterializeWorkflowDescriptorActor.{MaterializeWorkflowDescriptorCommand, MaterializeWorkflowDescriptorFailureResponse, MaterializeWorkflowDescriptorSuccessResponse}
import cromwell.util.SampleWdl.HelloWorld
import org.scalatest.BeforeAndAfter
import org.scalatest.mock.MockitoSugar
import spray.json.DefaultJsonProtocol._
import spray.json._
import wdl4s.values.{WdlInteger, WdlString}

import scala.concurrent.duration._
import scala.language.postfixOps

class MaterializeWorkflowDescriptorActorSpec extends CromwellTestkitSpec with BeforeAndAfter with MockitoSugar {

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

  val validInputsJson = HelloWorld.rawInputs.toJson.toString()
  val wdlSourceWithDocker = HelloWorld.wdlSource(""" runtime { docker: "ubuntu:latest" } """)
  val wdlSourceNoDocker = HelloWorld.wdlSource(""" runtime { } """)
  val Timeout = 10.second.dilated
  val NoBehaviourActor = system.actorOf(Props.empty)

  before {
  }

  after {
    system.stop(NoBehaviourActor)
  }

  "MaterializeWorkflowDescriptorActor" should {
    "accept valid WDL, inputs and options files" in {
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviourActor, workflowId))
      val sources = WorkflowSourceFiles(wdlSourceNoDocker, validInputsJson, validOptionsFile)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorSuccessResponse(wfDesc) =>
            wfDesc.id shouldBe workflowId
            wfDesc.name shouldBe "hello"
            wfDesc.namespace.tasks.size shouldBe 1
            wfDesc.workflowInputs.head shouldBe ("hello.hello.addressee", WdlString("world"))
            wfDesc.backendDescriptor.inputs.head shouldBe ("hello.hello.addressee", WdlString("world"))
            wfDesc.getWorkflowOption(WorkflowOptions.WriteToCache) shouldBe Some("true")
            wfDesc.getWorkflowOption(WorkflowOptions.ReadFromCache) shouldBe None
            // Default backend assignment is "Local":
            wfDesc.backendAssignments foreach {
              case (call, assignment) if call.task.name.equals("hello") => assignment shouldBe "Local"
              case (call, assignment) => fail(s"Unexpected call: ${call.task.name}")
            }
            wfDesc.engineFilesystems.size shouldBe 1
          case MaterializeWorkflowDescriptorFailureResponse(reason) => fail(s"Materialization failed with $reason")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    // Note to whoever comes next: I don't really know why this distinction exists. I've added this test but would
    // not be at all upset if the whole thing gets removed.
    "differently construct engine workflow inputs and backend inputs" in {
      val wdl =
        """
          |task bar { command { echo foobar } }
          |workflow foo {
          |  Int i
          |  Int j = 5
          |}
        """.stripMargin
      val inputs =
        """
          |{ "foo.i": "17" }
        """.stripMargin

      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviourActor, workflowId))
      val sources = WorkflowSourceFiles(wdl, inputs, validOptionsFile)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorSuccessResponse(wfDesc) =>


            wfDesc.workflowInputs foreach {
              case ("foo.i", wdlValue) => wdlValue shouldBe WdlInteger(17)
              case ("foo.j", wdlValue) => fail("Workflow declarations should not appear as workflow inputs")
              case (x, y) => fail(s"Unexpected input $x -> $y")
            }

            wfDesc.backendDescriptor.inputs foreach {
              case ("foo.i", wdlValue) => wdlValue shouldBe WdlInteger(17)
              case ("foo.j", wdlValue) => wdlValue shouldBe WdlInteger(5)
              case (x, y) => fail(s"Unexpected input $x -> $y")
            }
          case MaterializeWorkflowDescriptorFailureResponse(reason) => fail(s"Unexpected materialization failure: $reason")
          case unknown => fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }


    // TODO PBE: this should be done by MWDA (ticket #1076)
    "assign default runtime attributes" taggedAs PostMVP ignore {
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
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviourActor, workflowId))
      val sources = WorkflowSourceFiles(wdl, "{}", defaultDocker)
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
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviourActor, workflowId, cromwellBackends))
      val sources = WorkflowSourceFiles(wdl, "{}", "{}")
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, differentDefaultBackendConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorSuccessResponse(wfDesc) =>
            wfDesc.namespace.workflow.calls foreach {
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

      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviourActor, workflowId))
      val sources = WorkflowSourceFiles(wdl, "{}", "{}")
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, differentDefaultBackendConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            if (!reason.getMessage.contains("Backend for call foo.a ('NoSuchBackend') not registered in configuration file"))
              fail(s"Unexpected failure message from MaterializeWorkflowDescriptorActor: ${reason.getMessage}")
          case MaterializeWorkflowDescriptorSuccessResponse(wfDesc) => fail("This materialization should not have succeeded!")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "reject an invalid WDL source" in {
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviourActor, workflowId))
      val sources = WorkflowSourceFiles(unstructuredFile, validInputsJson, validOptionsFile)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed.\nUnable to load namespace from workflow: ERROR: Finished parsing without consuming all tokens.")
          case MaterializeWorkflowDescriptorSuccessResponse(wfDesc) => fail("This materialization should not have succeeded!")
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
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviourActor, workflowId))
      val sources = WorkflowSourceFiles(noWorkflowWdl, validInputsJson, validOptionsFile)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed.\nUnable to load namespace from workflow: Namespace does not have a local workflow to run")
          case MaterializeWorkflowDescriptorSuccessResponse(wfDesc) => fail("This materialization should not have succeeded!")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    // TODO: PBE: Re-enable (ticket #1063)
    "reject a taskless WDL source" taggedAs PostMVP ignore {
      val noWorkflowWdl =
        """
          |# no task foo { ... } block!!
          |
          |workflow foo {  }
        """.stripMargin
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviourActor, workflowId))
      val badWdlSources = WorkflowSourceFiles(noWorkflowWdl, validInputsJson, validOptionsFile)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(badWdlSources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed.\nUnable to load namespace from workflow: Namespace does not have a local workflow to run")
          case MaterializeWorkflowDescriptorSuccessResponse(wfDesc) => fail("This materialization should not have succeeded!")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }


    "reject an invalid options file" in {
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviourActor, workflowId))
      val sources = WorkflowSourceFiles(wdlSourceNoDocker, validInputsJson, unstructuredFile)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed.\nWorkflow contains invalid options JSON")
          case MaterializeWorkflowDescriptorSuccessResponse(wfDesc) => fail("This materialization should not have succeeded!")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "reject an invalid workflow inputs file" in {
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviourActor, workflowId))
      val sources = WorkflowSourceFiles(wdlSourceNoDocker, unstructuredFile, validOptionsFile)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed.\nWorkflow contains invalid inputs JSON")
          case MaterializeWorkflowDescriptorSuccessResponse(wfDesc) => fail("This materialization should not have succeeded!")
          case unknown =>
            fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }

    "reject requests if any required inputs are missing" in {
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviourActor, workflowId))
      val noInputsJson = "{}"
      val badOptionsSources = WorkflowSourceFiles(wdlSourceNoDocker, noInputsJson, validOptionsFile)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(badOptionsSources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed.\nRequired workflow input 'hello.hello.addressee' not specified")
          case MaterializeWorkflowDescriptorSuccessResponse(wfDesc) => fail("This materialization should not have succeeded!")
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
      val materializeWfActor = system.actorOf(MaterializeWorkflowDescriptorActor.props(NoBehaviourActor, workflowId))
      val sources = WorkflowSourceFiles(wdl, "{}", validOptionsFile)
      materializeWfActor ! MaterializeWorkflowDescriptorCommand(sources, minimumConf)

      within(Timeout) {
        expectMsgPF() {
          case MaterializeWorkflowDescriptorFailureResponse(reason) =>
            reason.getMessage should startWith("Workflow input processing failed.\nInvalid right-side type of 'foo.j'.  Expecting Int, got String")
          case MaterializeWorkflowDescriptorSuccessResponse(wfDesc) => fail("This materialization should not have succeeded!")
          case unknown => fail(s"Unexpected materialization response: $unknown")
        }
      }

      system.stop(materializeWfActor)
    }
  }
}
