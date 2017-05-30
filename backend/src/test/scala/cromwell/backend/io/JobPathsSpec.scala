package cromwell.backend.io

import com.typesafe.config.ConfigFactory
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptorKey, BackendSpec, TestConfig}
import cromwell.core.path.DefaultPathBuilder
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.TaskCall

class JobPathsSpec extends FlatSpec with Matchers with BackendSpec {

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

  val backendConfig =  ConfigFactory.parseString(configString)
  val defaultBackendConfigDescriptor = BackendConfigurationDescriptor(backendConfig, TestConfig.globalConfig)

  "JobPaths" should "provide correct paths for a job" in {

    val wd = buildWorkflowDescriptor(TestWorkflows.HelloWorld)
    val call: TaskCall = wd.workflow.taskCalls.head
    val jobKey = BackendJobDescriptorKey(call, None, 1)
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
    jobPaths.toDockerPath(DefaultPathBuilder.get(
      s"local-cromwell-executions/wf_hello/$id/call-hello/execution/stdout")).pathAsString shouldBe
      fullPath(s"/cromwell-executions/wf_hello/$id/call-hello/execution/stdout")
    jobPaths.toDockerPath(DefaultPathBuilder.get("/cromwell-executions/dock/path")).pathAsString shouldBe
      fullPath("/cromwell-executions/dock/path")

    val jobKeySharded = BackendJobDescriptorKey(call, Option(0), 1)
    val jobPathsSharded = new JobPathsWithDocker(workflowPaths, jobKeySharded)
    jobPathsSharded.callExecutionRoot.pathAsString shouldBe
      fullPath(s"local-cromwell-executions/wf_hello/$id/call-hello/shard-0/execution")

    val jobKeyAttempt = BackendJobDescriptorKey(call, None, 2)
    val jobPathsAttempt = new JobPathsWithDocker(workflowPaths, jobKeyAttempt)
    jobPathsAttempt.callExecutionRoot.pathAsString shouldBe
      fullPath(s"local-cromwell-executions/wf_hello/$id/call-hello/attempt-2/execution")

    val jobKeyShardedAttempt = BackendJobDescriptorKey(call, Option(0), 2)
    val jobPathsShardedAttempt = new JobPathsWithDocker(workflowPaths, jobKeyShardedAttempt)
    jobPathsShardedAttempt.callExecutionRoot.pathAsString shouldBe
      fullPath(s"local-cromwell-executions/wf_hello/$id/call-hello/shard-0/attempt-2/execution")
  }

  private def fullPath(path: String): String = DefaultPathBuilder.get(path).toAbsolutePath.pathAsString
}
