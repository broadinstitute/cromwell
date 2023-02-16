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

  it should "create the correct resources when an identity is passed in WorkflowOptions" in {
    val wei = Option("abc123")
    TesTask.makeResources(runtimeAttributes, wei) shouldEqual
        Resources(None, None, None, Option(false), None, Option(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> Option("abc123")))
    )
  }

  it should "create the correct resources when an empty identity is passed in WorkflowOptions" in {
    val wei = Option("")
    TesTask.makeResources(runtimeAttributes, wei) shouldEqual
        Resources(None, None, None, Option(false), None, Option(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> Option("")))
    )
  }

  it should "create the correct resources when no identity is passed in WorkflowOptions" in {
    val wei = None
    TesTask.makeResources(runtimeAttributes, wei) shouldEqual
        Resources(None, None, None, Option(false), None, Option(Map.empty[String, Option[String]])
    )
  }

  it should "create the correct resources when an identity is passed in via backend config" in {
    val weic = Option(WorkflowExecutionIdentityConfig("abc123"))
    val weio = Option(WorkflowExecutionIdentityOption("def456"))
    val wei = TesTask.getPreferredWorkflowExecutionIdentity(weic, weio)
    TesTask.makeResources(runtimeAttributes, wei) shouldEqual
        Resources(None, None, None, Option(false), None, Option(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> Option("abc123")))
    )
  }

  it should "create the correct resources when no identity is passed in via backend config" in {
    val weic = None
    val weio = Option(WorkflowExecutionIdentityOption("def456"))
    val wei = TesTask.getPreferredWorkflowExecutionIdentity(weic, weio)
    TesTask.makeResources(runtimeAttributes, wei) shouldEqual
        Resources(None, None, None, Option(false), None, Option(Map(TesWorkflowOptionKeys.WorkflowExecutionIdentity -> Option("def456")))
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

  it should "apply the correct transformation to a blob path" in {
    val orig = "https://coaexternalstorage.blob.core.windows.net/cromwell/mydir/myfile.txt"
    TesTask.transformBlobString(orig) shouldBe "/coaexternalstorage/cromwell/mydir/myfile.txt"
  }

  it should "not transform a non-blob path" in {
    val orig = "https://some-bogus-url.test/cromwell/mydir/myfile.txt"
    TesTask.transformBlobString(orig) shouldBe orig
  }

  it should "transform inputs" in {
    val baseInput = Input(
      Option("name"),
      Option("descr"),
      Option("https://coaexternalstorage.blob.core.windows.net/cromwell/mydir/myfile.txt"),
      "path",
      Option("type"),
      Option("content"),
    )
    val inputs = Option(Seq(baseInput))
    val outcome = inputs.map(TesTask.transformInputs)
    outcome.get.head.url.get shouldBe "/coaexternalstorage/cromwell/mydir/myfile.txt"
  }
}
