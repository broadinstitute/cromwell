package cromwell.backend.impl.tes

import common.assertion.CromwellTimeoutSpec
import common.mock.MockSugar
import cromwell.backend.validation.ContinueOnReturnCodeSet
import cromwell.backend.{BackendSpec, TestConfig}
import cromwell.core.WorkflowOptions
import cromwell.core.labels.Labels
import cromwell.core.logging.JobLogger
import cromwell.core.path.DefaultPathBuilder
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json.{JsObject, JsValue}
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
  val internalPathPrefix = ("internal_path_prefix", Option("mock/path/to/tes/task"))
  val additionalBackendParams = Map(internalPathPrefix)

  it should "create the correct resources when an identity is passed in WorkflowOptions" in {
    val wei = Option("abc123")
    TesTask.makeResources(runtimeAttributes, wei, additionalBackendParams) shouldEqual
        Resources(None, None, None, Option(false), None,
          Option(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> Option("abc123"),
            internalPathPrefix))
    )
  }

  it should "create the correct resources when an empty identity is passed in WorkflowOptions" in {
    val wei = Option("")
    TesTask.makeResources(runtimeAttributes, wei, additionalBackendParams) shouldEqual
      Resources(None, None, None, Option(false), None,
        Option(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> Option(""),
          internalPathPrefix))
    )
  }

  it should "create the correct resources when no identity is passed in WorkflowOptions" in {
    val wei = None
    TesTask.makeResources(runtimeAttributes, wei, additionalBackendParams) shouldEqual
      Resources(None, None, None, Option(false), None, Option(Map(internalPathPrefix)))
  }

  it should "create the correct resources when an identity is passed in via backend config" in {
    val weic = Option(WorkflowExecutionIdentityConfig("abc123"))
    val weio = Option(WorkflowExecutionIdentityOption("def456"))
    val wei = TesTask.getPreferredWorkflowExecutionIdentity(weic, weio)
    TesTask.makeResources(runtimeAttributes, wei, additionalBackendParams) shouldEqual
        Resources(None, None, None, Option(false), None, Option(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> Option("abc123"),
          internalPathPrefix))
    )
  }

  it should "create the correct resources when no identity is passed in via backend config" in {
    val weic = None
    val weio = Option(WorkflowExecutionIdentityOption("def456"))
    val wei = TesTask.getPreferredWorkflowExecutionIdentity(weic, weio)
    TesTask.makeResources(runtimeAttributes, wei, additionalBackendParams) shouldEqual
        Resources(None, None, None, Option(false), None, Option(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> Option("def456"),
          internalPathPrefix))
    )
  }

  it should "correctly set the internal path prefix when provided as a backend parameter" in {
    val wei = Option("abc123")
    val internalPathPrefix = ("internal_path_prefix", Option("mock/path/to/tes/task"))
    val additionalBackendParams = Map(internalPathPrefix)
    TesTask.makeResources(runtimeAttributes, wei, additionalBackendParams) shouldEqual
      Resources(None, None, None, Option(false), None,
        Option(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> Option("abc123"),
          internalPathPrefix))
      )
  }

  it should "correctly resolve the path to .../tes_task and add the k/v pair to backend parameters" in {
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

    //Assert path is created correctly
    val expectedKey = "internal_path_prefix"
    val expectedValue = Option(tesPaths.tesTaskRoot.pathAsString)

    //Assert path correctly ends up in the resources
    val additionalBackendParams = Map(expectedKey -> expectedValue)
    val wei = Option("abc123")
    TesTask.makeResources(runtimeAttributes, wei, additionalBackendParams) shouldEqual
      Resources(None, None, None, Option(false), None,
        Option(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> Option("abc123"),
          expectedKey -> expectedValue))
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
