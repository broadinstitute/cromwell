package cromwell.backend.impl.local

import cromwell.backend.BackendWorkflowDescriptor
import cromwell.backend.validation.ContinueOnReturnCodeSet
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.core.{WorkflowId, WorkflowOptions}
import org.scalatest.{Matchers, WordSpecLike}
import spray.json._
import wdl4s.WdlExpression.ScopedLookupFunction
import wdl4s.expression.NoFunctions
import wdl4s.util.TryUtil
import wdl4s.values.WdlValue
import wdl4s.{Call, NamespaceWithWorkflow, WdlExpression, WdlSource}

class LocalRuntimeAttributesSpec extends WordSpecLike with Matchers {

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
  val staticDefaults = new LocalRuntimeAttributes(ContinueOnReturnCodeSet(Set(0)), None, false)

  def workflowOptionsWithDefaultRA(defaults: Map[String, JsValue]) = {
    WorkflowOptions(JsObject(Map(
      "default_runtime_attributes" -> JsObject(defaults)
    )))
  }

  "LocalRuntimeAttributes" should {
    "return an instance of itself when there are no runtime attributes defined." in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      assertLocalRuntimeAttributesSuccessfulCreation(runtimeAttributes, emptyWorkflowOptions, staticDefaults)
    }

    "return an instance of itself when tries to validate a valid Docker entry" in {
      val expectedRuntimeAttributes = staticDefaults.copy(dockerImage = Option("ubuntu:latest"))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" }""").head
      assertLocalRuntimeAttributesSuccessfulCreation(runtimeAttributes, emptyWorkflowOptions, expectedRuntimeAttributes)
    }

    "return an instance of itself when tries to validate a valid Docker entry based on input" in {
      val expectedRuntimeAttributes = staticDefaults.copy(dockerImage = Option("you"))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "\${addressee}" }""").head
      assertLocalRuntimeAttributesSuccessfulCreation(runtimeAttributes, emptyWorkflowOptions, expectedRuntimeAttributes)
    }

    "use workflow options as default if docker key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      val workflowOptions = workflowOptionsWithDefaultRA(Map(DockerKey -> JsString("ubuntu:latest")))
      assertLocalRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions, staticDefaults.copy(dockerImage = Some("ubuntu:latest")))
    }

    "throw an exception when tries to validate an invalid Docker entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: 1 }""").head
      assertLocalRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting docker runtime attribute to be a String")
    }

    "return an instance of itself when tries to validate a valid failOnStderr entry" in {
      val expectedRuntimeAttributes = staticDefaults.copy(failOnStderr = true)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { failOnStderr: "true" }""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map(FailOnStderrKey -> JsBoolean(false)))
      assertLocalRuntimeAttributesSuccessfulCreation(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid failOnStderr entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { failOnStderr: "yes" }""").head
      assertLocalRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'")
    }

    "use workflow options as default if failOnStdErr key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      val workflowOptions = workflowOptionsWithDefaultRA(Map(FailOnStderrKey -> JsBoolean(true)))
      assertLocalRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions, staticDefaults.copy(failOnStderr = true))
    }

    "return an instance of itself when tries to validate a valid continueOnReturnCode entry" in {
      val expectedRuntimeAttributes = staticDefaults.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1)))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { continueOnReturnCode: 1 }""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map(ContinueOnReturnCodeKey -> JsBoolean(false)))
      assertLocalRuntimeAttributesSuccessfulCreation(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid continueOnReturnCode entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { continueOnReturnCode: "value" }""").head
      assertLocalRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]")
    }

    "use workflow options as default if continueOnReturnCode key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      val workflowOptions = workflowOptionsWithDefaultRA(Map(ContinueOnReturnCodeKey -> JsArray(Vector(JsNumber(1), JsNumber(2)))))
      assertLocalRuntimeAttributesSuccessfulCreation(runtimeAttributes, workflowOptions, staticDefaults.copy(continueOnReturnCode = ContinueOnReturnCodeSet(Set(1, 2))))
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
        val ra = call.task.runtimeAttributes.attrs mapValues { _.evaluate(createLookup(call), NoFunctions) }
        TryUtil.sequenceMap(ra, "Runtime attributes evaluation").get
    }
  }

  private def assertLocalRuntimeAttributesSuccessfulCreation(runtimeAttributes: Map[String, WdlValue], workflowOptions: WorkflowOptions, expectedRuntimeAttributes: LocalRuntimeAttributes): Unit = {
    try {
      assert(LocalRuntimeAttributes(runtimeAttributes, workflowOptions) == expectedRuntimeAttributes)
    } catch {
      case ex: RuntimeException => fail(s"Exception was not expected but received: ${ex.getMessage}")
    }
  }

  private def assertLocalRuntimeAttributesFailedCreation(runtimeAttributes: Map[String, WdlValue], exMsg: String): Unit = {
    try {
      LocalRuntimeAttributes(runtimeAttributes, WorkflowOptions(JsObject(Map.empty[String, JsValue])))
      fail("A RuntimeException was expected.")
    } catch {
      case ex: RuntimeException => assert(ex.getMessage.contains(exMsg))
    }
  }
}
