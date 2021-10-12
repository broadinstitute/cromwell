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
    val optionsMap: Map[String, JsValue] = excIdentity.map(i => TesWorkflowOptionKeys.Identity -> JsString(i)).toMap
    buildWdlWorkflowDescriptor(TestWorkflows.HelloWorld, None, WorkflowOptions(JsObject(optionsMap)))
  }

  "TesTask" should "create the correct resources when an identity is passed in WorkflowOptions" in {
    val wd = workflowDescriptorWithIdentity(Some("abc123"))
    assert(
      TesTask.makeResources(runtimeAttributes, wd)
        == Resources(None, None, None, Some(false), None, Some(Map(TesWorkflowOptionKeys.Identity -> "abc123")))
    )
  }

  "TesTask" should "create the correct resources when an empty identity is passed in WorkflowOptions" in {
    val wd = workflowDescriptorWithIdentity(Some(""))
    assert(
      TesTask.makeResources(runtimeAttributes, wd)
        == Resources(None, None, None, Some(false), None, Some(Map(TesWorkflowOptionKeys.Identity -> "")))
    )
  }

  "TesTask" should "create the correct resources when no identity is passed in WorkflowOptions" in {
    val wd = workflowDescriptorWithIdentity(None)
    assert(
      TesTask.makeResources(runtimeAttributes, wd)
        == Resources(None, None, None, Some(false), None, None)
    )
  }

  "TesTask" should "create the correct response when a numeric identity is passed in WorkflowOptions" in {
    val wd = buildWdlWorkflowDescriptor(
      TestWorkflows.HelloWorld,
      None,
      WorkflowOptions(
        JsObject(Map(TesWorkflowOptionKeys.Identity -> JsNumber(5)))
      )
    )
    assert(
      TesTask.makeResources(runtimeAttributes, wd)
        == Resources(None, None, None, Some(false), None, Some(Map(TesWorkflowOptionKeys.Identity -> "5")))
    )
  }
}
