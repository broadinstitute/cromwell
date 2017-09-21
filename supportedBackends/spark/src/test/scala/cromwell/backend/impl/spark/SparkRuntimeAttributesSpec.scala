package cromwell.backend.impl.spark

import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.{BackendWorkflowDescriptor, MemorySize}
import cromwell.core.labels.Labels
import cromwell.core.{NoIoFunctionSet, WorkflowId, WorkflowOptions}
import lenthall.validation.ErrorOr._
import org.scalatest.{Matchers, WordSpecLike}
import spray.json.{JsBoolean, JsNumber, JsObject, JsValue}
import wdl._
import wdl.values.WdlValue
import wom.executable.Executable.ResolvedExecutableInputs

class SparkRuntimeAttributesSpec extends WordSpecLike with Matchers {

  val HelloWorld =
    """
      |task hello {
      |  command {
      |    helloApp
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

  val staticDefaults = SparkRuntimeAttributes(1, MemorySize.parse("1 GB").get, None, "com.test.spark" , failOnStderr = false)

  def workflowOptionsWithDefaultRA(defaults: Map[String, JsValue]) = {
    WorkflowOptions(JsObject(Map(
      "default_runtime_attributes" -> JsObject(defaults)
    )))
  }

  "SparkRuntimeAttributes" should {
    "return an instance of itself when there are no runtime attributes defined." in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { %s: "%s" }""").head
      assertSparkRuntimeAttributes(runtimeAttributes, emptyWorkflowOptions, staticDefaults)
    }

    "return an instance of itself when tries to validate a valid Number of Executors entry" in {
      val expectedRuntimeAttributes = staticDefaults.copy(numberOfExecutors = Option(1))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { numberOfExecutors: 1 %s: "%s" }""").head
      assertSparkRuntimeAttributes(runtimeAttributes, emptyWorkflowOptions, expectedRuntimeAttributes)
    }

    "use workflow options as default if numberOfExecutors key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { %s: "%s" }""").head
      val workflowOptions = workflowOptionsWithDefaultRA(Map(SparkRuntimeAttributes.NumberOfExecutorsKey -> JsNumber(1)))
      assertSparkRuntimeAttributes(runtimeAttributes, workflowOptions, staticDefaults.copy(numberOfExecutors = Some(1)))
    }

    "throw an exception when tries to validate an invalid numberOfExecutors entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { numberOfExecutors: "1" %s: "%s" }""").head
      assertSparkRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting numberOfExecutors runtime attribute to be an Integer")
    }

    "return an instance of itself when tries to validate a valid failOnStderr entry" in {
      val expectedRuntimeAttributes = staticDefaults.copy(failOnStderr = true)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { failOnStderr: "true" %s: "%s" }""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map(FailOnStderrKey -> JsBoolean(false)))
      assertSparkRuntimeAttributes(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid failOnStderr entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { failOnStderr: "yes" %s: "%s" }""").head
      assertSparkRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'")
    }

    "use workflow options as default if failOnStdErr key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { %s: "%s" }""").head
      val workflowOptions = workflowOptionsWithDefaultRA(Map(FailOnStderrKey -> JsBoolean(true)))
      assertSparkRuntimeAttributes(runtimeAttributes, workflowOptions, staticDefaults.copy(failOnStderr = true))
    }

  }

  private def buildWorkflowDescriptor(wdl: WorkflowSource,
                                      inputs: ResolvedExecutableInputs = Map.empty,
                                      options: WorkflowOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue])),
                                      runtime: String) = {
    BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      WdlNamespaceWithWorkflow.load(wdl.replaceAll("RUNTIME", runtime.format("appMainClass", "com.test.spark")), Seq.empty[ImportResolver])
        .get.workflow.womDefinition.getOrElse(fail("Cannot build Wom Workflow")),
      inputs,
      options,
      Labels.empty
    )
  }

  private def createRuntimeAttributes(workflowSource: WorkflowSource, runtimeAttributes: String): List[Map[String, WdlValue]] = {  
    val workflowDescriptor = buildWorkflowDescriptor(workflowSource, runtime = runtimeAttributes)

    workflowDescriptor.workflow.taskCallNodes.toList map {
      call =>
        val staticValues = workflowDescriptor.knownValues.flatMap {
          case (outputPort, resolvedInput) => resolvedInput.select[WdlValue] map { outputPort.name -> _ }
        }
        val ra = call.callable.runtimeAttributes.attributes mapValues { _.evaluateValue(staticValues, NoIoFunctionSet) }
        ra.sequence.getOrElse(fail("Failed to evaluate runtime attributes"))
    }
  }

  private def assertSparkRuntimeAttributes(runtimeAttributes: Map[String, WdlValue], workflowOptions: WorkflowOptions, expectedRuntimeAttributes: SparkRuntimeAttributes): Unit = {
    try {
      assert(SparkRuntimeAttributes(runtimeAttributes, workflowOptions) == expectedRuntimeAttributes)
    } catch {
      case ex: RuntimeException => fail(s"Exception was not expected but received: ${ex.getMessage}")
    }
    ()
  }

  private def assertSparkRuntimeAttributesFailedCreation(runtimeAttributes: Map[String, WdlValue], exMsg: String): Unit = {
    try {
      SparkRuntimeAttributes(runtimeAttributes, emptyWorkflowOptions)
      fail("A RuntimeException was expected.")
    } catch {
      case ex: RuntimeException => assert(ex.getMessage.contains(exMsg))
    }
    ()
  }
}
