package cromwell.backend.impl.tes

import common.assertion.CromwellTimeoutSpec
import common.mock.MockSugar
import cromwell.backend.validation.ContinueOnReturnCodeSet
import cromwell.backend.{BackendJobDescriptor, BackendSpec, BackendWorkflowDescriptor, TestConfig}
import cromwell.core.{RootWorkflowId, WorkflowId, WorkflowOptions}
import cromwell.core.labels.Labels
import cromwell.core.logging.JobLogger
import cromwell.core.path.DefaultPathBuilder
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json.{JsObject, JsValue}
import wom.InstantiatedCommand

import java.util.UUID

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
  val internalPathPrefix = Option("mock/path/to/tes/task")
  val expectedTuple = "internal_path_prefix" -> internalPathPrefix

  it should "create the correct resources when an identity is passed in WorkflowOptions" in {
    val wei = Option("abc123")
    TesTask.makeResources(runtimeAttributes, wei, internalPathPrefix) shouldEqual
        Resources(None, None, None, Option(false), None,
          Option(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> Option("abc123"),
            expectedTuple))
    )
  }

  it should "create the correct resources when an empty identity is passed in WorkflowOptions" in {
    val wei = Option("")
    TesTask.makeResources(runtimeAttributes, wei, internalPathPrefix) shouldEqual
      Resources(None, None, None, Option(false), None,
        Option(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> Option(""),
          expectedTuple))
    )
  }

  it should "create the correct resources when no identity is passed in WorkflowOptions" in {
    val wei = None
    TesTask.makeResources(runtimeAttributes, wei, internalPathPrefix) shouldEqual
      Resources(None, None, None, Option(false), None, Option(Map(expectedTuple)))
  }

  it should "create the correct resources when an identity is passed in via backend config" in {
    val weic = Option(WorkflowExecutionIdentityConfig("abc123"))
    val weio = Option(WorkflowExecutionIdentityOption("def456"))
    val wei = TesTask.getPreferredWorkflowExecutionIdentity(weic, weio)
    TesTask.makeResources(runtimeAttributes, wei, internalPathPrefix) shouldEqual
        Resources(None, None, None, Option(false), None, Option(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> Option("abc123"),
          expectedTuple))
    )
  }

  it should "create the correct resources when no identity is passed in via backend config" in {
    val weic = None
    val weio = Option(WorkflowExecutionIdentityOption("def456"))
    val wei = TesTask.getPreferredWorkflowExecutionIdentity(weic, weio)
    TesTask.makeResources(runtimeAttributes, wei, internalPathPrefix) shouldEqual
        Resources(None, None, None, Option(false), None, Option(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> Option("def456"),
          expectedTuple))
    )
  }

  it should "correctly set the internal path prefix when provided as a backend parameter" in {
    val wei = Option("abc123")
    val internalPathPrefix = Option("mock/path/to/tes/task")
    TesTask.makeResources(runtimeAttributes, wei, internalPathPrefix) shouldEqual
      Resources(None, None, None, Option(false), None,
        Option(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> Option("abc123"),
          "internal_path_prefix" -> internalPathPrefix)
      ))
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

    val expectedKey = "internal_path_prefix"
    val expectedValue = Option(tesPaths.tesTaskRoot)

    //Assert path correctly ends up in the resources
    val wei = Option("abc123")
    TesTask.makeResources(runtimeAttributes, wei, expectedValue) shouldEqual
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

    task.tags shouldBe Option(
      Map(
        "foo" -> Option("bar"),
        "workflow_id" -> Option(workflowDescriptor.id.toString),
        "root_workflow_id" -> Option(workflowDescriptor.id.toString),
        "parent_workflow_id" -> None
      )
    )
  }

  it should "put workflow ids in tags" in {
    val jobLogger = mock[JobLogger]
    val emptyWorkflowOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue]))
    val workflowDescriptor = buildWdlWorkflowDescriptor(TestWorkflows.HelloWorld)
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

    task.tags shouldBe Option(
      Map(
        "workflow_id" -> Option(workflowDescriptor.id.toString),
        "root_workflow_id" -> Option(workflowDescriptor.id.toString),
        "parent_workflow_id" -> None
      )
    )
  }

  it should "put non-root workflow ids in tags" in {
    // Doing this test with mocks rather than real job/workflow descriptors as above because
    // getting the subworkflow structure build was really hard.
    val jobLogger = mock[JobLogger]
    val rootWorkflowId = RootWorkflowId(UUID.randomUUID())
    val subWorkflowId = WorkflowId(UUID.randomUUID())
    val subSubWorkflowId = WorkflowId(UUID.randomUUID())

    val jobDescriptor = mock[BackendJobDescriptor]
    val workflowDescriptor = mock[BackendWorkflowDescriptor]
    val tesPaths = mock[TesJobPaths]

    jobDescriptor.workflowDescriptor returns workflowDescriptor
    workflowDescriptor.customLabels returns Labels.empty
    workflowDescriptor.id returns subSubWorkflowId
    workflowDescriptor.rootWorkflowId returns rootWorkflowId
    workflowDescriptor.possibleParentWorkflowId returns Option(subWorkflowId)

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

    tesTask.tags shouldBe Map(
      "workflow_id" -> Option(subSubWorkflowId.toString),
      "root_workflow_id" -> Option(rootWorkflowId.toString),
      "parent_workflow_id" -> Option(subWorkflowId.toString)
    )
  }
}
