package cromwell.backend.impl.tes

import common.assertion.CromwellTimeoutSpec
import cromwell.backend.BackendSpec
import cromwell.backend.validation.ContinueOnReturnCodeSet
import cromwell.core.WorkflowOptions
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json.{JsObject, JsString, JsValue}

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
}
