package cromwell.backend.impl.sge

import cromwell.backend.BackendWorkflowDescriptor
import cromwell.backend.validation.RuntimeAttributesKeys._
import cromwell.backend.validation.{TryUtils, ContinueOnReturnCode, ContinueOnReturnCodeSet}
import cromwell.core.{WorkflowId, WorkflowOptions}
import org.scalatest.{Matchers, WordSpecLike}
import spray.json.{JsValue, JsObject}
import wdl4s.WdlExpression.ScopedLookupFunction
import wdl4s.expression.NoFunctions
import wdl4s.values.WdlValue
import wdl4s.{Call, WdlExpression, WdlSource, NamespaceWithWorkflow}

class SgeRuntimeAttributesSpec extends WordSpecLike with Matchers {

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
    Docker -> None,
    FailOnStderr -> false,
    ContinueOnReturnCode -> ContinueOnReturnCodeSet(Set(0)))

  "SgeRuntimeAttributes" should {
    "return an instance of itself when there are no runtime attributes defined." in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { }""").head
      assertSgeRuntimeAttributesSuccessfulCreation(runtimeAttributes, defaultRuntimeAttributes)
    }

    "return an instance of itself when tries to validate a valid Docker entry" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (Docker -> Option("ubuntu:latest"))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "ubuntu:latest" }""").head
      assertSgeRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "return an instance of itself when tries to validate a valid Docker entry based on input" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (Docker -> Option("you"))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: "\${addressee}" }""").head
      assertSgeRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid Docker entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { docker: 1 }""").head
      assertSgeRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting docker runtime attribute to be a String")
    }

    "return an instance of itself when tries to validate a valid failOnStderr entry" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (FailOnStderr -> true)
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { failOnStderr: "true" }""").head
      assertSgeRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid failOnStderr entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { failOnStderr: "yes" }""").head
      assertSgeRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting failOnStderr runtime attribute to be a Boolean or a String with values of 'true' or 'false'")
    }

    "return an instance of itself when tries to validate a valid continueOnReturnCode entry" in {
      val expectedRuntimeAttributes = defaultRuntimeAttributes + (ContinueOnReturnCode -> ContinueOnReturnCodeSet(Set(1)))
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { continueOnReturnCode: 1 }""").head
      assertSgeRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
    }

    "throw an exception when tries to validate an invalid continueOnReturnCode entry" in {
      val runtimeAttributes = createRuntimeAttributes(HelloWorld, """runtime { continueOnReturnCode: "value" }""").head
      assertSgeRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting continueOnReturnCode runtime attribute to be either a Boolean, a String 'true' or 'false', or an Array[Int]")
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
        TryUtils.sequenceMap(ra, "Runtime attributes evaluation").get
    }
  }

  private def assertSgeRuntimeAttributesSuccessfulCreation(runtimeAttributes: Map[String, WdlValue], expectedRuntimeAttributes: Map[String, Any]): Unit = {
    try {
      val sgeRuntimeAttributes = SgeRuntimeAttributes(runtimeAttributes)
      assert(sgeRuntimeAttributes.dockerImage == expectedRuntimeAttributes.get(Docker).get.asInstanceOf[Option[String]])
      assert(sgeRuntimeAttributes.failOnStderr == expectedRuntimeAttributes.get(FailOnStderr).get.asInstanceOf[Boolean])
      assert(sgeRuntimeAttributes.continueOnReturnCode == expectedRuntimeAttributes.get(ContinueOnReturnCode).get.asInstanceOf[ContinueOnReturnCode])
    } catch {
      case ex: RuntimeException => fail(s"Exception was not expected but received: ${ex.getMessage}")
    }
  }

  private def assertSgeRuntimeAttributesFailedCreation(runtimeAttributes: Map[String, WdlValue], exMsg: String): Unit = {
    try {
      SgeRuntimeAttributes(runtimeAttributes)
      fail("A RuntimeException was expected.")
    } catch {
      case ex: RuntimeException => assert(ex.getMessage.contains(exMsg))
    }
  }
}
