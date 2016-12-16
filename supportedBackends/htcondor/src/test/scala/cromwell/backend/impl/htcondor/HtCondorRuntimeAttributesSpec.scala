package cromwell.backend.impl.htcondor

import cromwell.backend.{BackendSpec, MemorySize}
import cromwell.backend.validation.ContinueOnReturnCodeSet
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.core.WorkflowOptions
import lenthall.util.TryUtil
import org.scalatest.{Matchers, WordSpecLike}
import spray.json._
import wdl4s.WdlExpression._
import wdl4s._
import wdl4s.expression.NoFunctions
import wdl4s.values.WdlValue

class HtCondorRuntimeAttributesSpec extends WordSpecLike with Matchers {

  import BackendSpec._

  val HelloWorld =
    s"""
      |task hello {
      |  String addressee = "you"
      |  command {
      |    echo "Hello $${addressee}!"
      |  }
      |  output {
      |    String salutation = read_string(stdout())
      |  }
      |
      |  RUNTIME
      |}
      |
      |workflow wf_hello {
      |  call hello
      |}
    """.stripMargin

  val emptyWorkflowOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue]))

  val memorySize = MemorySize.parse("0.512 GB").get
  val diskSize = MemorySize.parse("1.024 GB").get
  val staticDefaults = new HtCondorRuntimeAttributes(ContinueOnReturnCodeSet(Set(0)), None, None, None, false, 1,
    memorySize, diskSize, None)

  def workflowOptionsWithDefaultRA(defaults: Map[String, JsValue]) = {
    WorkflowOptions(JsObject(Map(
      "default_runtime_attributes" -> JsObject(defaults)
    )))
  }

  "HtCondorRuntimeAttributes" should {
    "return an instance of itself when there are no runtime attributes defined." in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, emptyWorkflowOptions, staticDefaults)
    }

    "return an instance of itself when tries to validate a valid Docker entry" in {
      val expectedRuntimeAttributes = staticDefaults.copy(dockerImage = Option("ubuntu:latest"))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" }""").head
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, emptyWorkflowOptions, expectedRuntimeAttributes)
    }

    "return an instance of itself when tries to validate a valid Docker entry based on input" in {
      val expectedRuntimeAttributes = staticDefaults.copy(dockerImage = Option("you"))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, s"""runtime { docker: "\\$${addressee}" }""").head
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, emptyWorkflowOptions, expectedRuntimeAttributes)
    }

    "use workflow options as default if docker key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      val workflowOptions = workflowOptionsWithDefaultRA(Map(DockerKey -> JsString("ubuntu:latest")))
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions, staticDefaults.copy(dockerImage = Some("ubuntu:latest")))
    }

    "throw an exception when tries to validate an invalid Docker entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: 1 }""").head
      assertHtCondorRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting docker runtime attribute to be a String")
    }

    "return an instance of itself when tries to validate a valid docker working directory entry" in {
      val expectedRuntimeAttributes = staticDefaults.copy(dockerWorkingDir = Option("/workingDir"))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { dockerWorkingDir: "/workingDir" }""").head
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, emptyWorkflowOptions, expectedRuntimeAttributes)
    }

    "return an instance of itself when tries to validate a valid docker working directory entry based on input" in {
      val expectedRuntimeAttributes = staticDefaults.copy(dockerWorkingDir = Option("you"))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, s"""runtime { dockerWorkingDir: "\\$${addressee}" }""").head
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, emptyWorkflowOptions, expectedRuntimeAttributes)
    }

    "use workflow options as default if docker working directory key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      val workflowOptions = workflowOptionsWithDefaultRA(Map("dockerWorkingDir" -> JsString("/workingDir")))
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions, staticDefaults.copy(dockerWorkingDir = Some("/workingDir")))
    }

    "throw an exception when tries to validate an invalid docker working directory entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { dockerWorkingDir: 1 }""").head
      assertHtCondorRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting dockerWorkingDir runtime attribute to be a String")
    }

    "return an instance of itself when tries to validate a valid docker output directory entry" in {
      val expectedRuntimeAttributes = staticDefaults.copy(dockerOutputDir = Option("/outputDir"))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { dockerOutputDir: "/outputDir" }""").head
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, emptyWorkflowOptions, expectedRuntimeAttributes)
    }

    "return an instance of itself when tries to validate a valid docker output directory entry based on input" in {
      val expectedRuntimeAttributes = staticDefaults.copy(dockerOutputDir = Option("you"))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, s"""runtime { dockerOutputDir: "\\$${addressee}" }""").head
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, emptyWorkflowOptions, expectedRuntimeAttributes)
    }

    "use workflow options as default if docker output directory key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      val workflowOptions = workflowOptionsWithDefaultRA(Map("dockerOutputDir" -> JsString("/outputDir")))
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions, staticDefaults.copy(dockerOutputDir = Some("/outputDir")))
    }

    "throw an exception when tries to validate an invalid docker output directory entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { dockerOutputDir: 1 }""").head
      assertHtCondorRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting dockerOutputDir runtime attribute to be a String")
    }

    "return an instance of itself when tries to validate a valid failOnStderr entry" in {
      val expectedRuntimeAttributes = staticDefaults.copy(failOnStderr = true)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { failOnStderr: "true" }""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map(FailOnStderrKey -> JsBoolean(false)))
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid failOnStderr entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { failOnStderr: "yes" }""").head
      assertHtCondorRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'")
    }

    "use workflow options as default if failOnStdErr key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      val workflowOptions = workflowOptionsWithDefaultRA(Map(FailOnStderrKey -> JsBoolean(true)))
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions, staticDefaults.copy(failOnStderr = true))
    }

    "return an instance of itself when tries to validate a valid continueOnReturnCode entry" in {
      val expectedRuntimeAttributes = staticDefaults.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1)))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { continueOnReturnCode: 1 }""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map(ContinueOnReturnCodeKey -> JsBoolean(false)))
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid continueOnReturnCode entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { continueOnReturnCode: "value" }""").head
      assertHtCondorRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]")
    }

    "use workflow options as default if continueOnReturnCode key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      val workflowOptions = workflowOptionsWithDefaultRA(Map(ContinueOnReturnCodeKey -> JsArray(Vector(JsNumber(1), JsNumber(2)))))
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions, staticDefaults.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1, 2))))
    }

    "return an instance of itself when tries to validate a valid cpu entry" in {
      val expectedRuntimeAttributes = staticDefaults.copy(cpu = 2)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { cpu: 2 }""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map(CpuKey -> JsString("6")))
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid cpu entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { cpu: "value" }""").head
      assertHtCondorRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting cpu runtime attribute to be an Integer")
    }

    "use workflow options as default if cpu key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      val expectedRuntimeAttributes = staticDefaults.copy(cpu = 6)
      val workflowOptions = workflowOptionsWithDefaultRA(Map(CpuKey -> JsString("6")))
      val workflowOptions2 = workflowOptionsWithDefaultRA(Map(CpuKey -> JsNumber("6")))
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions, expectedRuntimeAttributes)
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions2, expectedRuntimeAttributes)
    }

    "use default cpu value when there is no cpu key entry" in {
      val expectedRuntimeAttributes = staticDefaults.copy(cpu = 1)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map())
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "return an instance of itself when tries to validate a valid memory entry" in {
      val expectedRuntimeAttributes = staticDefaults.copy(memory = MemorySize.parse("1 GB").get)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { memory: "1 GB" }""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map(MemoryKey -> JsString("blahaha")))
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid memory entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" memory: "value" }""").head
      assertHtCondorRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting memory runtime attribute to be an Integer or String with format '8 GB'")
    }

    "use workflow options as default if memory key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      val expectedRuntimeAttributes = staticDefaults.copy(memory = MemorySize.parse("65 GB").get)
      val workflowOptions = workflowOptionsWithDefaultRA(Map(MemoryKey -> JsString("65 GB")))
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions, expectedRuntimeAttributes)
    }

    "use default memory value when there is no memory key entry" in {
      val expectedRuntimeAttributes = staticDefaults.copy(memory = MemorySize.parse("0.512 GB").get)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map())
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "return an instance of itself when tries to validate a valid disk entry" in {
      val expectedRuntimeAttributes = staticDefaults.copy(disk = MemorySize.parse("1 GB").get)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { disk: "1 GB" }""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map("disk" -> JsString("blahaha")))
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid String disk entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" disk: "value" }""").head
      assertHtCondorRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting memory runtime attribute to be an Integer or String with format '8 GB'")
    }

    "throw an exception when tries to validate an invalid Integer array disk entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" disk: [1] }""").head
      assertHtCondorRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting disk runtime attribute to be an Integer or String with format '8 GB'")
    }

    "use workflow options as default if disk key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      val expectedRuntimeAttributes = staticDefaults.copy(disk = MemorySize.parse("65 GB").get)
      val workflowOptions = workflowOptionsWithDefaultRA(Map("disk" -> JsString("65 GB")))
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions, expectedRuntimeAttributes)
    }

    "use default disk value when there is no disk key entry" in {
      val expectedRuntimeAttributes = staticDefaults.copy(disk = MemorySize.parse("1.024 GB").get)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map())
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "return an instance of itself when tries to validate a valid native specs entry" in {
      val expectedRuntimeAttributes = staticDefaults.copy(nativeSpecs = Option(Array("spec1", "spec2")))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { nativeSpecs: ["spec1", "spec2"] }""").head
      assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes, emptyWorkflowOptions, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid native specs entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { nativeSpecs: [1, 2] }""").head
      assertHtCondorRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting nativeSpecs runtime attribute to be an Array of Strings.")
    }
  }

  private def assertHtCondorRuntimeAttributesSuccessfulCreation(runtimeAttributes: Map[String, WdlValue],
                                                                workflowOptions: WorkflowOptions,
                                                                expectedRuntimeAttributes: HtCondorRuntimeAttributes) = {
    try {
      val actualRuntimeAttr = HtCondorRuntimeAttributes(runtimeAttributes, workflowOptions)
      assert(actualRuntimeAttr.cpu == expectedRuntimeAttributes.cpu)
      assert(actualRuntimeAttr.disk == expectedRuntimeAttributes.disk)
      assert(actualRuntimeAttr.memory == expectedRuntimeAttributes.memory)
      assert(actualRuntimeAttr.continueOnReturnCode == expectedRuntimeAttributes.continueOnReturnCode)
      assert(actualRuntimeAttr.failOnStderr == expectedRuntimeAttributes.failOnStderr)
      assert(actualRuntimeAttr.dockerWorkingDir == expectedRuntimeAttributes.dockerWorkingDir)
      assert(actualRuntimeAttr.dockerImage == expectedRuntimeAttributes.dockerImage)
      assert(actualRuntimeAttr.dockerOutputDir == expectedRuntimeAttributes.dockerOutputDir)
      expectedRuntimeAttributes.nativeSpecs match {
        case Some(ns) => assert(ns.deep == expectedRuntimeAttributes.nativeSpecs.get.deep)
        case None => assert(expectedRuntimeAttributes.nativeSpecs.isEmpty)
      }
    } catch {
      case ex: RuntimeException => fail(s"Exception was not expected but received: ${ex.getMessage}")
    }
  }

  private def assertHtCondorRuntimeAttributesFailedCreation(runtimeAttributes: Map[String, WdlValue], exMsg: String) = {
    try {
      HtCondorRuntimeAttributes(runtimeAttributes, emptyWorkflowOptions)
      fail("A RuntimeException was expected.")
    } catch {
      case ex: RuntimeException => assert(ex.getMessage.contains(exMsg))
    }
  }

  private def createRuntimeAttributes(wdlSource: WdlSource, runtimeAttributes: String): Seq[Map[String, WdlValue]] = {
    val workflowDescriptor = buildWorkflowDescriptor(wdlSource, runtime = runtimeAttributes)

    def createLookup(call: Call): ScopedLookupFunction = {
      val knownInputs = workflowDescriptor.knownValues
      call.lookupFunction(knownInputs, NoFunctions)
    }

    workflowDescriptor.workflow.taskCalls.toSeq map {
      call =>
        val ra = call.task.runtimeAttributes.attrs mapValues { _.evaluate(createLookup(call), NoFunctions) }
        TryUtil.sequenceMap(ra, "Runtime attributes evaluation").get
    }
  }
}
