package cromwell.backend.impl.jes

import cromwell.backend.impl.jes.io.{DiskType, JesAttachedDisk, JesWorkingDisk}
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.validation.{ContinueOnReturnCodeFlag, ContinueOnReturnCodeSet}
import cromwell.backend.{BackendWorkflowDescriptor, MemorySize}
import cromwell.core.{WorkflowId, WorkflowOptions}
import org.scalatest.{Matchers, WordSpecLike}
import spray.json._
import wdl4s.WdlExpression.ScopedLookupFunction
import wdl4s.expression.NoFunctions
import wdl4s.parser.MemoryUnit
import wdl4s.util.TryUtil
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

  val emptyWorkflowOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue]))

  def workflowOptionsWithDefaultRA(defaults: Map[String, JsValue]) = {
    WorkflowOptions(JsObject(Map(
      "default_runtime_options" -> JsObject(defaults)
    )))
  }

  val staticDefaults = new JesRuntimeAttributes(1, Vector("us-central1-a"), 0, 10, MemorySize(2, MemoryUnit.GB), Seq(JesWorkingDisk(DiskType.SSD, 10)), None, false, ContinueOnReturnCodeSet(Set(0)))
  val staticDefaultsWithUbuntu = staticDefaults.copy(dockerImage = Some("ubuntu:latest"))

  "JesRuntimeAttributes" should {

    "throw an exception when there are no runtime attributes defined." in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Can't find an attribute value for key docker")
    }

    "use workflow options as default if docker key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      val workflowOptions = workflowOptionsWithDefaultRA(Map(DockerKey -> JsString("ubuntu:latest")))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions, staticDefaults.copy(dockerImage = Some("ubuntu:latest")))
    }

    "return an instance of itself when tries to validate a valid Docker entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" }""").head
      val expectedRuntimeAttributes = staticDefaults.copy(dockerImage = Option("ubuntu:latest"))
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map(DockerKey -> JsString("ubuntu:fromoptions")))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "return an instance of itself when tries to validate a valid Docker entry based on input" in {
      val expectedRuntimeAttributes = staticDefaults.copy(dockerImage = Option("you"))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "\${addressee}" }""").head
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, emptyWorkflowOptions, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid Docker entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: 1 }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting docker runtime attribute to be a String")
    }

    "return an instance of itself when tries to validate a valid failOnStderr entry" in {
      val expectedRuntimeAttributes = staticDefaultsWithUbuntu.copy(failOnStderr = true)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" failOnStderr: "true" }""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map(FailOnStderrKey -> JsString("false")))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid failOnStderr entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" failOnStderr: "yes" }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'")
    }

    "use workflow options as default if failOnStderr key is missing" in {
      val expectedRuntimeAttributes = staticDefaultsWithUbuntu.copy(failOnStderr = true)
      val workflowOptions1 = workflowOptionsWithDefaultRA(Map(FailOnStderrKey -> JsString("true")))
      val workflowOptions2 = workflowOptionsWithDefaultRA(Map(FailOnStderrKey -> JsBoolean(true)))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" failOnStderr: "true" }""").head
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions1, expectedRuntimeAttributes)
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions2, expectedRuntimeAttributes)
    }

    "return an instance of itself when tries to validate a valid continueOnReturnCode entry" in {
      val expectedRuntimeAttributes = staticDefaultsWithUbuntu.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1)))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" continueOnReturnCode: 1 }""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map(ContinueOnReturnCodeKey -> JsBoolean(true)))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid continueOnReturnCode entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" continueOnReturnCode: "value" }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]")
    }

    "use workflow options as default if continueOnReturnCode key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest"}""").head

      val expectedRuntimeAttributesBool = staticDefaultsWithUbuntu.copy(continueOnReturnCode = ContinueOnReturnCodeFlag(true))
      val workflowOptionsBool1 = workflowOptionsWithDefaultRA(Map(ContinueOnReturnCodeKey -> JsString("true")))
      val workflowOptionsBool2 = workflowOptionsWithDefaultRA(Map(ContinueOnReturnCodeKey -> JsBoolean(true)))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptionsBool1, expectedRuntimeAttributesBool)
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptionsBool2, expectedRuntimeAttributesBool)

      val expectedRuntimeAttributesSet = staticDefaultsWithUbuntu.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1, 2, 3)))
      val workflowOptionsSet = workflowOptionsWithDefaultRA(Map(ContinueOnReturnCodeKey -> JsArray(Vector(JsNumber(1), JsNumber(2), JsNumber(3)))))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptionsSet, expectedRuntimeAttributesSet)
    }

    "return an instance of itself when tries to validate a valid cpu entry" in {
      val expectedRuntimeAttributes = staticDefaultsWithUbuntu.copy(cpu = 2)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" cpu: 2 }""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map(CpuKey -> JsString("6")))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid cpu entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" cpu: "value" }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting cpu runtime attribute to be an Integer")
    }

    "use workflow options as default if cpu key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest"}""").head

      val expectedRuntimeAttributes = staticDefaultsWithUbuntu.copy(cpu = 6)
      val workflowOptions = workflowOptionsWithDefaultRA(Map(CpuKey -> JsString("6")))
      val workflowOptions2 = workflowOptionsWithDefaultRA(Map(CpuKey -> JsNumber("6")))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions, expectedRuntimeAttributes)
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions2, expectedRuntimeAttributes)
    }

    "return an instance of itself when tries to validate a valid zones entry" in {
      val expectedRuntimeAttributes = staticDefaultsWithUbuntu.copy(zones = Vector("us-central1-z"))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" zones: "us-central1-z" }""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map(JesRuntimeAttributes.ZonesKey -> JsString("zone-from-options")))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid zones entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" zones: 1 }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting zones runtime attribute to be either a whitespace separated String or an Array[String]")
    }

    "return InitializationSuccess when tries to validate a valid array zones entry" in {
      val expectedRuntimeAttributes = staticDefaultsWithUbuntu.copy(zones = Vector("us-central1-y", "us-central1-z"))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" zones: ["us-central1-y", "us-central1-z"] }""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map(JesRuntimeAttributes.ZonesKey -> JsString("zone-from-options")))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid array zones entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" zones: [2, 1] }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting zones runtime attribute to be either a whitespace separated String or an Array[String]")
    }

    "use workflow options as default if zones key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest"}""").head

      val expectedRuntimeAttributes = staticDefaultsWithUbuntu.copy(zones = Vector("us-central1-y"))
      val workflowOptions = workflowOptionsWithDefaultRA(Map(JesRuntimeAttributes.ZonesKey -> JsString("us-central1-y")))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions, expectedRuntimeAttributes)

      val expectedRuntimeAttributes2 = staticDefaultsWithUbuntu.copy(zones = Vector("us-central1-y", "us-central1-z"))
      val workflowOptions2 = workflowOptionsWithDefaultRA(Map(JesRuntimeAttributes.ZonesKey -> JsArray(Vector(JsString("us-central1-y"), JsString("us-central1-z")))))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions2, expectedRuntimeAttributes2)
    }

    "return an instance of itself when tries to validate a valid preemptible entry" in {
      val expectedRuntimeAttributes = staticDefaultsWithUbuntu.copy(preemptible = 3)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" preemptible: 3 }""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map(JesRuntimeAttributes.PreemptibleKey -> JsString("4")))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid preemptible entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" preemptible: "value" }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting preemptible runtime attribute to be an Integer")
    }

    "use workflow options as default if preemptible key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest"}""").head

      val expectedRuntimeAttributes = staticDefaultsWithUbuntu.copy(preemptible = 4)
      val workflowOptions = workflowOptionsWithDefaultRA(Map(JesRuntimeAttributes.PreemptibleKey -> JsString("4")))
      val workflowOptions2 = workflowOptionsWithDefaultRA(Map(JesRuntimeAttributes.PreemptibleKey -> JsNumber("4")))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions, expectedRuntimeAttributes)
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions2, expectedRuntimeAttributes)
    }

    "return InitializationSuccess when tries to validate a valid bootDiskSizeGb entry" in {
      val expectedRuntimeAttributes = staticDefaultsWithUbuntu.copy(bootDiskSize = 4)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" bootDiskSizeGb: 4 }""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map(JesRuntimeAttributes.BootDiskSizeKey -> JsString("10")))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid bootDiskSizeGb entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" bootDiskSizeGb: "value" }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting bootDiskSizeGb runtime attribute to be an Integer")
    }

    "use workflow options as default if bootDiskSizeGb key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest"}""").head

      val expectedRuntimeAttributes = staticDefaultsWithUbuntu.copy(bootDiskSize = 50)
      val workflowOptions = workflowOptionsWithDefaultRA(Map(JesRuntimeAttributes.BootDiskSizeKey -> JsString("50")))
      val workflowOptions2 = workflowOptionsWithDefaultRA(Map(JesRuntimeAttributes.BootDiskSizeKey -> JsNumber("50")))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions, expectedRuntimeAttributes)
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions2, expectedRuntimeAttributes)
    }

    "return an instance of itself when tries to validate a valid disks entry" in {
      val expectedRuntimeAttributes = staticDefaultsWithUbuntu.copy(disks = Seq(JesAttachedDisk.parse("local-disk 20 SSD").get))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" disks: "local-disk 20 SSD" }""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map(JesRuntimeAttributes.DisksKey -> JsString("blahaha")))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid disks entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" disks: 10 }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting disks runtime attribute to be a comma separated String or Array[String]")
    }

    "throw an exception when tries to validate an invalid array disks entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" disks: [10, 11] }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting disks runtime attribute to be a comma separated String or Array[String]")
    }

    "use workflow options as default if disks key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest"}""").head

      val expectedRuntimeAttributes = staticDefaultsWithUbuntu.copy(disks = Seq(JesAttachedDisk.parse("local-disk 50 SSD").get))
      val workflowOptions = workflowOptionsWithDefaultRA(Map(JesRuntimeAttributes.DisksKey -> JsString("local-disk 50 SSD")))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions, expectedRuntimeAttributes)

      val expectedRuntimeAttributes2 = staticDefaultsWithUbuntu.copy(disks = Seq(JesAttachedDisk.parse("local-disk 50 SSD").get, JesAttachedDisk.parse("/mount/point 80 SSD").get))
      val workflowOptions2 = workflowOptionsWithDefaultRA(Map(JesRuntimeAttributes.DisksKey -> JsArray(Vector(JsString("local-disk 50 SSD"), JsString("/mount/point 80 SSD")))))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions2, expectedRuntimeAttributes2)
    }

    "return an instance of itself when tries to validate a valid memory entry" in {
      val expectedRuntimeAttributes = staticDefaultsWithUbuntu.copy(memory = MemorySize.parse("1 GB").get)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" memory: "1 GB" }""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map(MemoryKey -> JsString("blahaha")))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid memory entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" memory: "value" }""").head
      assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting memory runtime attribute to be an Integer or String with format '8 GB'")
    }

    "use workflow options as default if memory key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest"}""").head

      val expectedRuntimeAttributes = staticDefaultsWithUbuntu.copy(memory = MemorySize.parse("65 GB").get)
      val workflowOptions = workflowOptionsWithDefaultRA(Map(MemoryKey -> JsString("65 GB")))
      assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions, expectedRuntimeAttributes)
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
        TryUtil.sequenceMap(ra, "Runtime attributes evaluation").get
    }
  }

  private def assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes: Map[String, WdlValue], options: WorkflowOptions, expectedRuntimeAttributes: JesRuntimeAttributes): Unit = {
    try {
      assert(JesRuntimeAttributes(runtimeAttributes, options) == expectedRuntimeAttributes)
    } catch {
      case ex: RuntimeException => fail(s"Exception was not expected but received: ${ex.getMessage}")
    }
  }

  private def assertJesRuntimeAttributesFailedCreation(runtimeAttributes: Map[String, WdlValue], exMsg: String): Unit = {
    try {
      JesRuntimeAttributes(runtimeAttributes, emptyWorkflowOptions)
      fail("A RuntimeException was expected.")
    } catch {
      case ex: RuntimeException => assert(ex.getMessage.contains(exMsg))
    }
  }
}
