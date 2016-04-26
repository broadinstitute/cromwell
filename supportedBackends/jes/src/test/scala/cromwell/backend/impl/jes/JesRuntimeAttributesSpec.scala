package cromwell.backend.impl.jes

import cromwell.backend.BackendWorkflowDescriptor
import cromwell.backend.impl.jes.io.JesAttachedDisk
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.validation.{ContinueOnReturnCode, ContinueOnReturnCodeSet, MemorySize, TryUtils}
import cromwell.core.{WorkflowId, WorkflowOptions}
import org.scalatest.{Matchers, WordSpecLike}
import spray.json.{JsObject, JsValue}
import wdl4s.WdlExpression.ScopedLookupFunction
import wdl4s.expression.NoFunctions
import wdl4s.parser.MemoryUnit
import wdl4s.values.WdlValue
import wdl4s.{Call, NamespaceWithWorkflow, WdlExpression, WdlSource}

class JesRuntimeAttributesSpec extends WordSpecLike with Matchers {

  val HelloWorld =
    """
      |task hello {
      |  String addressee = "you"
      |  command {
      |    echo "Hello ${addressee}!"
      |  }
      |  output {
      |    String salutation = read_string(stdout())
      |  }
      |
      |  RUNTIME
      |}
      |
      |workflow hello {
      |  call hello
      |}
    """.stripMargin

  val defaultRuntimeAttributes = Map(
    Docker -> Some("ubuntu:latest"),
    FailOnStderr -> false,
    ContinueOnReturnCode -> ContinueOnReturnCodeSet(Set(0)),
    Cpu -> 1,
    Memory -> MemorySize(2.0, MemoryUnit.GB),
    JesRuntimeAttributes.ZonesKey -> JesRuntimeAttributes.ZoneDefaultValue.toVector,
    JesRuntimeAttributes.BootDiskSizeKey -> JesRuntimeAttributes.BootDiskSizeDefaultValue,
    JesRuntimeAttributes.PreemptibleKey -> JesRuntimeAttributes.PreemptibleDefaultValue,
    JesRuntimeAttributes.DisksKey -> Seq(JesRuntimeAttributes.DefaultJesWorkingDisk)
  )

  "JesRuntimeAttributes" should {

    "throw an exception when there are no runtime attributes defined." in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Failed to get Docker mandatory key from runtime attributes")
    }

    "return an instance of itself when tries to validate a valid Docker entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" }""").head
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, defaultRuntimeAttributes)
    }

    "return an instance of itself when tries to validate a valid Docker entry based on input" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (Docker -> Option("you"))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "\${addressee}" }""").head
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid Docker entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: 1 }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting docker runtime attribute to be a String")
    }

    "return an instance of itself when tries to validate a valid failOnStderr entry" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (FailOnStderr -> true)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" failOnStderr: "true" }""").head
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid failOnStderr entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" failOnStderr: "yes" }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'")
    }

    "return an instance of itself when tries to validate a valid continueOnReturnCode entry" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (ContinueOnReturnCode -> ContinueOnReturnCodeSet(Set(1)))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" failOnStderr: "false" continueOnReturnCode: 1 }""").head
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid continueOnReturnCode entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" failOnStderr: "yes" continueOnReturnCode: "value" }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]")
    }

    "return an instance of itself when tries to validate a valid cpu entry" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (Cpu -> 2)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" cpu: 2 }""").head
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid cpu entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" cpu: "value" }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting cpu runtime attribute to be an Integer")
    }

    "return an instance of itself when tries to validate a valid zones entry" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (JesRuntimeAttributes.ZonesKey -> Seq("us-central1-z").toVector)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" zones: "us-central1-z" }""").head
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid zones entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" zones: 1 }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting zones runtime attribute to be either a whitespace separated String or an Array[String]")
    }

    "return InitializationSuccess when tries to validate a valid array zones entry" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (JesRuntimeAttributes.ZonesKey -> Seq("us-central1-y", "us-central1-z").toVector)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" zones: ["us-central1-y", "us-central1-z"] }""").head
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid array zones entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" zones: [2, 1] }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting zones runtime attribute to be either a whitespace separated String or an Array[String]")
    }

    "return an instance of itself when tries to validate a valid preemptible entry" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (JesRuntimeAttributes.PreemptibleKey -> 3)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" preemptible: 3 }""").head
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid preemptible entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" preemptible: "value" }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting preemptible runtime attribute to be an Integer")
    }

    "return InitializationSuccess when tries to validate a valid bootDiskSizeGb entry" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (JesRuntimeAttributes.BootDiskSizeKey -> 4)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" bootDiskSizeGb: 4 }""").head
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid bootDiskSizeGb entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" bootDiskSizeGb: "value" }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting bootDiskSizeGb runtime attribute to be an Integer")
    }

    "return an instance of itself when tries to validate a valid disks entry" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (JesRuntimeAttributes.DisksKey -> Seq(JesAttachedDisk.parse("local-disk 20 SSD").get))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" disks: "local-disk 20 SSD" }""").head
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid disks entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" disks: 10 }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting disks runtime attribute to be a comma separated String or Array[String]")
    }

    "throw an exception when tries to validate an invalid array disks entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" disks: [10, 11] }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting disks runtime attribute to be a comma separated String or Array[String]")
    }

    "return an instance of itself when tries to validate a valid memory entry" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (Memory -> MemorySize.parse("1 GB").get)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" memory: "1 GB" }""").head
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid memory entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" memory: "value" }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting memory runtime attribute to be an Integer or String with format '8 GB'")
    }
  }

  private def buildWorkflowDescriptor(wdl: WdlSource,
                                      inputs: Map[String, WdlValue] = Map.empty,
                                      options: WorkflowOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue])),
                                      runtime: String = "") = {
    new BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      NamespaceWithWorkflow.load(wdl.replaceAll("RUNTIME", runtime)),
      inputs,
      options
    )
  }

  private def createRuntimeAttributes(wdlSource: WdlSource, runtimeAttributes: String = "") = {
    val workflowDescriptor = buildWorkflowDescriptor(wdlSource, runtime = runtimeAttributes)

    def createLookup(call: Call): ScopedLookupFunction = {
      val declarations = workflowDescriptor.workflowNamespace.workflow.declarations ++ call.task.declarations
      val knownInputs = workflowDescriptor.inputs
      WdlExpression.standardLookupFunction(knownInputs, declarations, NoFunctions)
    }

    workflowDescriptor.workflowNamespace.workflow.calls map {
      call =>
        val ra = call.task.runtimeAttributes.attrs mapValues {
          _.evaluate(createLookup(call), NoFunctions)
        }
        TryUtils.sequenceMap(ra, "Runtime attributes evaluation").get
    }
  }

  private def assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes: Map[String, WdlValue], expectedRuntimeAttributes: Map[String, Any]): Unit = {
    try {
      val jesRuntimeAttributes = JesRuntimeAttributes(runtimeAttributes)
      assert(jesRuntimeAttributes.dockerImage == expectedRuntimeAttributes.get(Docker).get.asInstanceOf[Option[String]])
      assert(jesRuntimeAttributes.failOnStderr == expectedRuntimeAttributes.get(FailOnStderr).get.asInstanceOf[Boolean])
      assert(jesRuntimeAttributes.continueOnReturnCode == expectedRuntimeAttributes.get(ContinueOnReturnCode).get.asInstanceOf[ContinueOnReturnCode])
      assert(jesRuntimeAttributes.cpu == expectedRuntimeAttributes.get(Cpu).get.asInstanceOf[Int])
      assert(jesRuntimeAttributes.memory == expectedRuntimeAttributes.get(Memory).get.asInstanceOf[MemorySize])
      assert(jesRuntimeAttributes.zones == expectedRuntimeAttributes.get(JesRuntimeAttributes.ZonesKey).get.asInstanceOf[Vector[String]])
      assert(jesRuntimeAttributes.bootDiskSize == expectedRuntimeAttributes.get(JesRuntimeAttributes.BootDiskSizeKey).get.asInstanceOf[Int])
      assert(jesRuntimeAttributes.preemptible == expectedRuntimeAttributes.get(JesRuntimeAttributes.PreemptibleKey).get.asInstanceOf[Int])
      assert(jesRuntimeAttributes.disks == expectedRuntimeAttributes.get(JesRuntimeAttributes.DisksKey).get.asInstanceOf[Seq[JesAttachedDisk]])
    } catch {
      case ex: RuntimeException => fail(s"Exception was not expected but received: ${ex.getMessage}")
    }
  }

  private def assertJesRuntimeAttributesFailedCreation(runtimeAttributes: Map[String, WdlValue], exMsg: String): Unit = {
    try {
      JesRuntimeAttributes(runtimeAttributes)
      fail("A RuntimeException was expected.")
    } catch {
      case ex: RuntimeException => assert(ex.getMessage.contains(exMsg))
    }
  }
}
