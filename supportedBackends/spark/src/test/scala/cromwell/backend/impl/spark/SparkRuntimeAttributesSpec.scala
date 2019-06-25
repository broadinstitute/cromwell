package cromwell.backend.impl.spark

import common.collections.EnhancedCollections._
import common.validation.ErrorOr._
import cromwell.backend.BackendWorkflowDescriptor
import cromwell.core.labels.Labels
import cromwell.core.{HogGroup, WorkflowId, WorkflowOptions}
import org.scalatest.{Matchers, WordSpecLike}
import spray.json.{JsBoolean, JsNumber, JsObject, JsString, JsValue}
import wdl.draft2.model.{Draft2ImportResolver, WdlNamespaceWithWorkflow}
import wdl.transforms.draft2.wdlom2wom._
import wom.RuntimeAttributesKeys._
import wom.core.WorkflowSource
import wom.expression.NoIoFunctionSet
import wom.graph.GraphNodePort.OutputPort
import wom.transforms.WomWorkflowDefinitionMaker.ops._
import wom.values.WomValue
import wom.format.MemorySize

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

  val staticDefaults = SparkRuntimeAttributes(1, MemorySize.parse("1 GB").get, None, None, None, failOnStderr = false)

  def workflowOptionsWithDefaultRA(defaults: Map[String, JsValue]) = {
    WorkflowOptions(JsObject(Map(
      "default_runtime_attributes" -> JsObject(defaults)
    )))
  }

  "SparkRuntimeAttributes" should {
    "return an instance of itself when there are no runtime attributes defined." in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime {}""").head
      assertSparkRuntimeAttributes(runtimeAttributes, emptyWorkflowOptions, staticDefaults)
    }

    "return an instance of itself when tries to validate a valid Number of Executors entry" in {
      val expectedRuntimeAttributes = staticDefaults.copy(numberOfExecutors = Option(1))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { numberOfExecutors: 1}""").head
      assertSparkRuntimeAttributes(runtimeAttributes, emptyWorkflowOptions, expectedRuntimeAttributes)
    }

    "use workflow options as default if numberOfExecutors key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime {}""").head
      val workflowOptions = workflowOptionsWithDefaultRA(Map(SparkRuntimeAttributes.NumberOfExecutorsKey -> JsNumber(1)))
      assertSparkRuntimeAttributes(runtimeAttributes, workflowOptions, staticDefaults.copy(numberOfExecutors = Some(1)))
    }

    "throw an exception when tries to validate an invalid numberOfExecutors entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { numberOfExecutors: "1"}""").head
      assertSparkRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting numberOfExecutors runtime attribute to be an Integer")
    }

    "return an instance of itself when tries to validate a valid failOnStderr entry" in {
      val expectedRuntimeAttributes = staticDefaults.copy(failOnStderr = true)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { failOnStderr: "true"}""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map(FailOnStderrKey -> JsBoolean(false)))
      assertSparkRuntimeAttributes(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid failOnStderr entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { failOnStderr: "yes"}""").head
      assertSparkRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'")
    }

    "use workflow options as default if failOnStdErr key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime {}""").head
      val workflowOptions = workflowOptionsWithDefaultRA(Map(FailOnStderrKey -> JsBoolean(true)))
      assertSparkRuntimeAttributes(runtimeAttributes, workflowOptions, staticDefaults.copy(failOnStderr = true))
    }

    "return an instance of itself when tries to validate a valid appMainClass entry" in {
      val expectedRuntimeAttributes = staticDefaults.copy(appMainClass = Some("com.test.spark"))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { appMainClass: "com.test.spark"}""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map(SparkRuntimeAttributes.AppMainClassKey -> JsString("com.test.spark")))
      assertSparkRuntimeAttributes(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid appMainClass entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { appMainClass: 1}""").head
      assertSparkRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting appMainClass runtime attribute to be a String")
    }

    "use workflow options as default if appMainClass key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime {}""").head
      val workflowOptions = workflowOptionsWithDefaultRA(Map(SparkRuntimeAttributes.AppMainClassKey -> JsString("com.test.spark")))
      assertSparkRuntimeAttributes(runtimeAttributes, workflowOptions, staticDefaults.copy(appMainClass = Some("com.test.spark")))
    }

    "return an instance of itself when tries to validate a valid additionalArgs entry" in {
      val expectedRuntimeAttributes = staticDefaults.copy(additionalArgs = Some("--a 'a1' --b 'b1' --c 'c1 d1=-d2'"))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { additionalArgs: "--a 'a1' --b 'b1' --c 'c1 d1=-d2'"}""").head
      val shouldBeIgnored = workflowOptionsWithDefaultRA(Map(SparkRuntimeAttributes.AdditionalArgsKey -> JsString("--a 'a1' --b 'b1' --c 'c1 d1=-d2'")))
      assertSparkRuntimeAttributes(runtimeAttributes, shouldBeIgnored, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid additionalArgs entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { additionalArgs: 1}""").head
      assertSparkRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting additionalArgs runtime attribute to be a String")
    }

    "use workflow options as default if additionalArgs key is missing" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime {}""").head
      val workflowOptions = workflowOptionsWithDefaultRA(Map(SparkRuntimeAttributes.AdditionalArgsKey -> JsString("--a 'a1' --b 'b1' --c 'c1 d1=-d2'")))
      assertSparkRuntimeAttributes(runtimeAttributes, workflowOptions, staticDefaults.copy(additionalArgs = Some( "--a 'a1' --b 'b1' --c 'c1 d1=-d2'")))
    }

  }

  private def buildWorkflowDescriptor(wdl: WorkflowSource,
                                      inputs: Map[OutputPort, WomValue] = Map.empty,
                                      options: WorkflowOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue])),
                                      runtime: String) = {
    val wdlNamespace = WdlNamespaceWithWorkflow.load(wdl.replaceAll("RUNTIME", runtime), Seq.empty[Draft2ImportResolver]).get

    BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      wdlNamespace.workflow.toWomWorkflowDefinition(isASubworkflow = false).getOrElse(fail("Cannot build Wom Workflow")),
      inputs,
      options,
      Labels.empty,
      HogGroup("foo"),
      List.empty,
      None
    )
  }

  private def createRuntimeAttributes(workflowSource: WorkflowSource, runtimeAttributes: String): List[Map[String, WomValue]] = {
    val workflowDescriptor = buildWorkflowDescriptor(workflowSource, runtime = runtimeAttributes)

    workflowDescriptor.callable.taskCallNodes.toList map {
      call =>
        val staticValues = workflowDescriptor.knownValues.map {
          case (outputPort, resolvedInput) => outputPort.name -> resolvedInput
        }
        val ra = call.callable.runtimeAttributes.attributes safeMapValues { _.evaluateValue(staticValues, NoIoFunctionSet) }
        ra.sequence.getOrElse(fail("Failed to evaluate runtime attributes"))
    }
  }

  private def assertSparkRuntimeAttributes(runtimeAttributes: Map[String, WomValue], workflowOptions: WorkflowOptions, expectedRuntimeAttributes: SparkRuntimeAttributes): Unit = {
    try {
      assert(SparkRuntimeAttributes(runtimeAttributes, workflowOptions) == expectedRuntimeAttributes)
    } catch {
      case ex: RuntimeException => fail(s"Exception was not expected but received: ${ex.getMessage}")
    }
    ()
  }

  private def assertSparkRuntimeAttributesFailedCreation(runtimeAttributes: Map[String, WomValue], exMsg: String): Unit = {
    try {
      SparkRuntimeAttributes(runtimeAttributes, emptyWorkflowOptions)
      fail("A RuntimeException was expected.")
    } catch {
      case ex: RuntimeException => assert(ex.getMessage.contains(exMsg))
    }
    ()
  }
}
