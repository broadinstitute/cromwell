package cromwell.backend.io

import com.typesafe.config.ConfigFactory
import common.assertion.CromwellTimeoutSpec
import cromwell.backend.io.JobPathsSpecHelper._
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptorKey, BackendSpec, TestConfig}
import cromwell.core.path.DefaultPathBuilder
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wom.graph.CommandCallNode

class JobPathsSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with BackendSpec {

  val configString =
    """
      |        root: "local-cromwell-executions"
      |        filesystems {
      |          local {
      |            localization: [
      |              "hard-link", "soft-link", "copy"
      |            ]
      |          }
      |          gcs {
      |            auth = "application-default"
      |          }
      |        }
    """.stripMargin

  val backendConfig = ConfigFactory.parseString(configString)
  val defaultBackendConfigDescriptor = BackendConfigurationDescriptor(backendConfig, TestConfig.globalConfig)

  "JobPaths" should "provide correct paths for a job" in {

    val wd = buildWdlWorkflowDescriptor(TestWorkflows.HelloWorld)
    val call: CommandCallNode = wd.callable.taskCallNodes.head
    val jobKey: BackendJobDescriptorKey = BackendJobDescriptorKey(call, None, 1)
    val workflowPaths = new WorkflowPathsWithDocker(wd, backendConfig)
    val jobPaths = new JobPathsWithDocker(workflowPaths, jobKey)
    val id = wd.id
    jobPaths.callRoot.pathAsString shouldBe
      fullPath(s"local-cromwell-executions/wf_hello/$id/call-hello")
    jobPaths.callExecutionRoot.pathAsString shouldBe
      fullPath(s"local-cromwell-executions/wf_hello/$id/call-hello/execution")
    jobPaths.returnCode.pathAsString shouldBe
      fullPath(s"local-cromwell-executions/wf_hello/$id/call-hello/execution/rc")
    jobPaths.script.pathAsString shouldBe
      fullPath(s"local-cromwell-executions/wf_hello/$id/call-hello/execution/script")
    jobPaths.stderr.pathAsString shouldBe
      fullPath(s"local-cromwell-executions/wf_hello/$id/call-hello/execution/stderr")
    jobPaths.stdout.pathAsString shouldBe
      fullPath(s"local-cromwell-executions/wf_hello/$id/call-hello/execution/stdout")
    jobPaths.callExecutionRoot.pathAsString shouldBe
      fullPath(s"local-cromwell-executions/wf_hello/$id/call-hello/execution")
    jobPaths.callDockerRoot.pathAsString shouldBe
      fullPath(s"/cromwell-executions/wf_hello/$id/call-hello")
    jobPaths.callExecutionDockerRoot.pathAsString shouldBe
      fullPath(s"/cromwell-executions/wf_hello/$id/call-hello/execution")
    jobPaths
      .toDockerPath(DefaultPathBuilder.get(s"local-cromwell-executions/wf_hello/$id/call-hello/execution/stdout"))
      .pathAsString shouldBe
      fullPath(s"/cromwell-executions/wf_hello/$id/call-hello/execution/stdout")
    jobPaths.toDockerPath(DefaultPathBuilder.get("/cromwell-executions/dock/path")).pathAsString shouldBe
      fullPath("/cromwell-executions/dock/path")
    jobPaths.memoryRetryRC.pathAsString shouldBe
      fullPath(s"local-cromwell-executions/wf_hello/$id/call-hello/execution/memory_retry_rc")

    val jobKeySharded: BackendJobDescriptorKey = BackendJobDescriptorKey(call, Option(0), 1)
    val jobPathsSharded = new JobPathsWithDocker(workflowPaths, jobKeySharded)
    jobPathsSharded.callExecutionRoot.pathAsString shouldBe
      fullPath(s"local-cromwell-executions/wf_hello/$id/call-hello/shard-0/execution")

    val jobKeyAttempt: BackendJobDescriptorKey = BackendJobDescriptorKey(call, None, 2)
    val jobPathsAttempt = new JobPathsWithDocker(workflowPaths, jobKeyAttempt)
    jobPathsAttempt.callExecutionRoot.pathAsString shouldBe
      fullPath(s"local-cromwell-executions/wf_hello/$id/call-hello/attempt-2/execution")

    val jobKeyShardedAttempt: BackendJobDescriptorKey = BackendJobDescriptorKey(call, Option(0), 2)
    val jobPathsShardedAttempt = new JobPathsWithDocker(workflowPaths, jobKeyShardedAttempt)
    jobPathsShardedAttempt.callExecutionRoot.pathAsString shouldBe
      fullPath(s"local-cromwell-executions/wf_hello/$id/call-hello/shard-0/attempt-2/execution")
  }

  private def fullPath(path: String): String = DefaultPathBuilder.get(path).toAbsolutePath.pathAsString
}
