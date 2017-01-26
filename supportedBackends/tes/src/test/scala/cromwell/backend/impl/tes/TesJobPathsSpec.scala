package cromwell.backend.impl.tes

import better.files._
import com.typesafe.config.ConfigFactory
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptorKey, BackendSpec}
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.TaskCall

class TesJobPathsSpec extends FlatSpec with Matchers with BackendSpec {

  val configString =
    """
      |        root = "local-cromwell-executions"
      |        dockerRoot = "/cromwell-executions"
      |
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

  val globalConfig = ConfigFactory.load()
  val backendConfig =  ConfigFactory.parseString(configString)
  val defaultBackendConfigDescriptor = BackendConfigurationDescriptor(backendConfig, globalConfig)

  "JobPaths" should "provide correct paths for a job" in {

    val wd = buildWorkflowDescriptor(TestWorkflows.HelloWorld)
    val call: TaskCall = wd.workflow.taskCalls.head
    val jobKey = BackendJobDescriptorKey(call, None, 1)
    val jobPaths = new TesJobPaths(jobKey, wd, backendConfig)
    val id = wd.id
    jobPaths.callRoot.toString shouldBe
      File(s"local-cromwell-executions/wf_hello/$id/call-hello").pathAsString
    jobPaths.callExecutionRoot.toString shouldBe
      File(s"local-cromwell-executions/wf_hello/$id/call-hello/execution").pathAsString
    jobPaths.returnCode.toString shouldBe
      File(s"local-cromwell-executions/wf_hello/$id/call-hello/execution/rc").pathAsString
    jobPaths.script.toString shouldBe
      File(s"local-cromwell-executions/wf_hello/$id/call-hello/execution/script").pathAsString
    jobPaths.stderr.toString shouldBe
      File(s"local-cromwell-executions/wf_hello/$id/call-hello/execution/stderr").pathAsString
    jobPaths.stdout.toString shouldBe
      File(s"local-cromwell-executions/wf_hello/$id/call-hello/execution/stdout").pathAsString
    jobPaths.callExecutionRoot.toString shouldBe
      File(s"local-cromwell-executions/wf_hello/$id/call-hello/execution").pathAsString
    jobPaths.callDockerRoot.toString shouldBe
      File(s"/cromwell-executions/wf_hello/$id/call-hello").pathAsString
    jobPaths.callExecutionDockerRoot.toString shouldBe
      File(s"/cromwell-executions/wf_hello/$id/call-hello/execution").pathAsString

    val jobKeySharded = BackendJobDescriptorKey(call, Option(0), 1)
    val jobPathsSharded = new TesJobPaths(jobKeySharded, wd, backendConfig)
    jobPathsSharded.callExecutionRoot.toString shouldBe
      File(s"local-cromwell-executions/wf_hello/$id/call-hello/shard-0/execution").pathAsString

    val jobKeyAttempt = BackendJobDescriptorKey(call, None, 2)
    val jobPathsAttempt = new TesJobPaths(jobKeyAttempt, wd, backendConfig)
    jobPathsAttempt.callExecutionRoot.toString shouldBe
      File(s"local-cromwell-executions/wf_hello/$id/call-hello/attempt-2/execution").pathAsString

    val jobKeyShardedAttempt = BackendJobDescriptorKey(call, Option(0), 2)
    val jobPathsShardedAttempt = new TesJobPaths(jobKeyShardedAttempt, wd, backendConfig)
    jobPathsShardedAttempt.callExecutionRoot.toString shouldBe
      File(s"local-cromwell-executions/wf_hello/$id/call-hello/shard-0/attempt-2/execution").pathAsString
  }
}
