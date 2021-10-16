package cromwell.backend.impl.tes

import common.assertion.CromwellTimeoutSpec
import cromwell.backend.BackendSpec
import cromwell.backend.validation.ContinueOnReturnCodeSet
import cromwell.core.WorkflowOptions
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json.{JsNumber, JsObject, JsString, JsValue}

class TesTaskSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with BackendSpec {

  val runtimeAttributes = new TesRuntimeAttributes(
    ContinueOnReturnCodeSet(Set(0)),
    "ubuntu:latest",
    None,
    false,
    None,
    None,
    None,
    false
  )

  def workflowDescriptorWithIdentity(excIdentity: Option[String]) = {
    val optionsMap: Map[String, JsValue] = excIdentity.map(i => TesWorkflowOptionKeys.WorkflowExecutionIdentity -> JsString(i)).toMap
    buildWdlWorkflowDescriptor(TestWorkflows.HelloWorld, None, WorkflowOptions(JsObject(optionsMap)))
  }

  it should "create the correct resources when an identity is passed in WorkflowOptions" in {
    val wd = workflowDescriptorWithIdentity(Option("abc123"))
    assert(
      TesTask.makeResources(runtimeAttributes, wd)
        == Resources(None, None, None, Option(false), None, Option(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> "abc123")))
    )
  }

  it should "create the correct resources when an empty identity is passed in WorkflowOptions" in {
    val wd = workflowDescriptorWithIdentity(Option(""))
    assert(
      TesTask.makeResources(runtimeAttributes, wd)
        == Resources(None, None, None, Option(false), None, Option(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> "")))
    )
  }

  it should "create the correct resources when no identity is passed in WorkflowOptions" in {
    val wd = workflowDescriptorWithIdentity(None)
    assert(
      TesTask.makeResources(runtimeAttributes, wd)
        == Resources(None, None, None, Option(false), None, Option(Map.empty[String, String]))
    )
  }

  it should "create the correct response when a numeric identity is passed in WorkflowOptions" in {
    val wd = buildWdlWorkflowDescriptor(
      TestWorkflows.HelloWorld,
      None,
      WorkflowOptions(
        JsObject(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> JsNumber(5)))
      )
    )
    assert(
      TesTask.makeResources(runtimeAttributes, wd)
        == Resources(None, None, None, Option(false), None, Option(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> "5")))
    )
  }

  // TODO this isn't actually the behavior we want
  it should "silently do nothing when the identity passed in WorkflowOptions is an object" in {
    val wd = buildWdlWorkflowDescriptor(
      TestWorkflows.HelloWorld,
      None,
      WorkflowOptions(
        JsObject(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> JsObject(Map("hi" -> JsString("there")))))
      )
    )
    assert(
      TesTask.makeResources(runtimeAttributes, wd)
        == Resources(None, None, None, Option(false), None, Option(Map.empty[String, String]))
    )
  }
}
