package cromwell.backend.io

import com.typesafe.config.Config
import cromwell.backend.{BackendJobBreadCrumb, BackendSpec, BackendWorkflowDescriptor}
import cromwell.core.path.DefaultPathBuilder
import cromwell.core.{JobKey, WorkflowId}
import cromwell.util.WomMocks
import org.mockito.Mockito._
import org.scalatest.{FlatSpec, Matchers}
import wom.graph.WomIdentifier

class WorkflowPathsSpec extends FlatSpec with Matchers with BackendSpec {

  def createConfig(values: Map[String, String]): Config = {
    val config = mock[Config]
    values.foreach {
      case (key: String, value: String) =>
        when(config.hasPath(key)).thenReturn(true)
        when(config.getString(key)).thenReturn(value)
    }
    config
  }

  def rootConfig(root: Option[String], dockerRoot: Option[String]) = {
    val values: Map[String,String] = root.map("root" -> _).toMap ++ dockerRoot.map("dockerRoot" -> _).toMap
    createConfig(values)
  }

  case class TestConfig(root: Option[String], dockerRoot: Option[String])
  val testConfigs: List[TestConfig] = List(
    TestConfig(None, None), // Defaults testing
    TestConfig(Some("local-cromwell-executions"), None), // Testing only with defined root
    TestConfig(None, Some("/dockerRootExecutions")), // Testing only with defined dockerRoot
    TestConfig(Some("local-cromwell-executions"), Some("/dockerRootExecutions")) // Testing with both defined.
  )

  def testWorkflowPaths(root: Option[String], dockerRoot: Option[String]) = {
    val backendConfig = rootConfig(root, dockerRoot)
    val wd = buildWdlWorkflowDescriptor(TestWorkflows.HelloWorld)
    val workflowPaths = new WorkflowPathsWithDocker(wd, backendConfig)
    val id = wd.id
    val expectedRoot = root.getOrElse(WorkflowPaths.DefaultExecutionRootString)
    val expectedDockerRoot = dockerRoot.getOrElse(WorkflowPathsWithDocker.DefaultDockerRoot)
    workflowPaths.workflowRoot.pathAsString shouldBe
      DefaultPathBuilder.get(s"$expectedRoot/wf_hello/$id").toAbsolutePath.pathAsString
    workflowPaths.dockerWorkflowRoot.pathAsString shouldBe
      s"$expectedDockerRoot/wf_hello/$id"
  }

  def testSubWorkflowPaths(root: Option[String], dockerRoot: Option[String]) = {
    val backendConfig = rootConfig(root, dockerRoot)
    val rootWd = mock[BackendWorkflowDescriptor]
    val rootWorkflow = WomMocks.mockWorkflowDefinition("rootWorkflow")
    val rootWorkflowId = WorkflowId.randomId()
    rootWd.callable returns rootWorkflow
    rootWd.id returns rootWorkflowId

    val subWd = mock[BackendWorkflowDescriptor]
    val subWorkflow = WomMocks.mockWorkflowDefinition("subWorkflow")
    val subWorkflowId = WorkflowId.randomId()
    subWd.callable returns subWorkflow
    subWd.id returns subWorkflowId

    val call1 = WomMocks.mockTaskCall(WomIdentifier("call1"))

    val jobKey = new JobKey {
      override def node = call1
      override def tag: String = "tag1"
      override def index: Option[Int] = Option(1)
      override def attempt: Int = 2
    }

    subWd.breadCrumbs returns List(BackendJobBreadCrumb(rootWorkflow, rootWorkflowId, jobKey))
    subWd.id returns subWorkflowId

    val workflowPaths = new WorkflowPathsWithDocker(subWd, backendConfig)
    val expectedRoot = root.getOrElse(WorkflowPaths.DefaultExecutionRootString)
    val expectedDockerRoot = dockerRoot.getOrElse(WorkflowPathsWithDocker.DefaultDockerRoot)

    workflowPaths.workflowRoot.pathAsString shouldBe
      DefaultPathBuilder.get(
        s"$expectedRoot/rootWorkflow/$rootWorkflowId/call-call1/shard-1/attempt-2/subWorkflow/$subWorkflowId"
      ).toAbsolutePath.pathAsString
    workflowPaths.dockerWorkflowRoot.pathAsString shouldBe s"$expectedDockerRoot/rootWorkflow/$rootWorkflowId/call-call1/shard-1/attempt-2/subWorkflow/$subWorkflowId"
  }

  "WorkflowPaths" should "provide correct paths for a workflow" in {
    testConfigs.foreach(config => testWorkflowPaths(config.root, config.dockerRoot))
  }

  "WorkflowPaths" should "provide correct paths for a sub workflow" in {
    testConfigs.foreach(config => testSubWorkflowPaths(config.root, config.dockerRoot))
  }
}
