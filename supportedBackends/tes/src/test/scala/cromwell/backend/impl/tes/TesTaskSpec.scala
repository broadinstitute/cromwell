package cromwell.backend.impl.tes

import common.assertion.CromwellTimeoutSpec
import common.mock.MockSugar
import cromwell.backend.{BackendSpec, TestConfig}
import cromwell.backend.validation.ContinueOnReturnCodeSet
import cromwell.core.WorkflowOptions
import cromwell.core.labels.Labels
import cromwell.core.logging.JobLogger
import cromwell.core.path.DefaultPathBuilder
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json.{JsObject, JsString, JsValue}
import wom.InstantiatedCommand

class TesTaskSpec
  extends AnyFlatSpec
    with CromwellTimeoutSpec
    with Matchers
    with BackendSpec
    with MockSugar {

  val runtimeAttributes = new TesRuntimeAttributes(
    ContinueOnReturnCodeSet(Set(0)),
    "ubuntu:latest",
    None,
    false,
    None,
    None,
    None,
    false,
    Map.empty
  )

  def workflowDescriptorWithIdentity(excIdentity: Option[String]) = {
    val optionsMap: Map[String, JsValue] = excIdentity.map(i => TesWorkflowOptionKeys.WorkflowExecutionIdentity -> JsString(i)).toMap
    buildWdlWorkflowDescriptor(TestWorkflows.HelloWorld, None, WorkflowOptions(JsObject(optionsMap)))
  }

  it should "create the correct resources when an identity is passed in WorkflowOptions" in {
    val wd = workflowDescriptorWithIdentity(Option("abc123"))
    TesTask.makeResources(runtimeAttributes, wd) shouldEqual
        Resources(None, None, None, Option(false), None, Option(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> Option("abc123")))
    )
  }

  it should "create the correct resources when an empty identity is passed in WorkflowOptions" in {
    val wd = workflowDescriptorWithIdentity(Option(""))
    TesTask.makeResources(runtimeAttributes, wd) shouldEqual
        Resources(None, None, None, Option(false), None, Option(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> Option("")))
    )
  }

  it should "create the correct resources when no identity is passed in WorkflowOptions" in {
    val wd = workflowDescriptorWithIdentity(None)
    TesTask.makeResources(runtimeAttributes, wd) shouldEqual
        Resources(None, None, None, Option(false), None, Option(Map.empty[String, Option[String]])
    )
  }


  it should "copy labels to tags" in {
    val jobLogger = mock[JobLogger]
    val emptyWorkflowOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue]))
    val workflowDescriptor = buildWdlWorkflowDescriptor(TestWorkflows.HelloWorld,
                                                        labels = Labels("foo" -> "bar"))
    val jobDescriptor = jobDescriptorFromSingleCallWorkflow(workflowDescriptor,
                                                            Map.empty,
                                                            emptyWorkflowOptions,
                                                            Set.empty)
    val tesPaths = TesJobPaths(jobDescriptor.key,
                               jobDescriptor.workflowDescriptor,
                               TestConfig.emptyConfig)
    val tesTask = TesTask(jobDescriptor,
                          TestConfig.emptyBackendConfigDescriptor,
                          jobLogger,
                          tesPaths,
                          runtimeAttributes,
                          DefaultPathBuilder.build("").get,
                          "",
                          InstantiatedCommand("command"),
                          "",
                          Map.empty,
                          "",
                          OutputMode.ROOT)

    val task = TesTask.makeTask(tesTask)

    task.tags shouldBe Option(Map("foo" -> "bar"))
  }
}
