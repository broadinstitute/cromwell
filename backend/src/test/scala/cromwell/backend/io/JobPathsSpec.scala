package cromwell.backend.io

import java.nio.file.Paths

import com.typesafe.config.ConfigFactory
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptorKey, BackendSpec}
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.Call

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

  val globalConfig = ConfigFactory.load()
  val backendConfig =  ConfigFactory.parseString(configString)
  val defaultBackendConfigDescriptor = new BackendConfigurationDescriptor(backendConfig, globalConfig)

  "JobPaths" should "provide correct paths for a job" in {

    val wd = buildWorkflowDescriptor(TestWorkflows.HelloWorld)
    val call: Call = wd.workflowNamespace.workflow.calls.head
    val jobKey = new BackendJobDescriptorKey(call, None, 1)
    val jobPaths = new JobPaths(wd, backendConfig, jobKey, None)
    val id = wd.id
    jobPaths.callRoot.toString shouldBe
      s"local-cromwell-executions/hello/$id/call-hello"
    jobPaths.returnCode.toString shouldBe
      s"local-cromwell-executions/hello/$id/call-hello/rc"
    jobPaths.script.toString shouldBe
      s"local-cromwell-executions/hello/$id/call-hello/script"
    jobPaths.stderr.toString shouldBe
      s"local-cromwell-executions/hello/$id/call-hello/stderr"
    jobPaths.stdout.toString shouldBe
      s"local-cromwell-executions/hello/$id/call-hello/stdout"
    jobPaths.callRoot.toString shouldBe
      s"local-cromwell-executions/hello/$id/call-hello"
    jobPaths.callDockerRoot.toString shouldBe
      s"/root/hello/$id/call-hello"
    jobPaths.toDockerPath(Paths.get(s"local-cromwell-executions/hello/$id/call-hello/stdout")).toString shouldBe
      s"/root/hello/$id/call-hello/stdout"
    jobPaths.toDockerPath(Paths.get("/root/dock/path")).toString shouldBe
      "/root/dock/path"

    val jobKeySharded = new BackendJobDescriptorKey(call, Option(0), 1)
    val jobPathsSharded = new JobPaths(wd, backendConfig, jobKeySharded, None)
    jobPathsSharded.callRoot.toString shouldBe
      s"local-cromwell-executions/hello/$id/call-hello/shard-0"

    val jobKeyAttempt = new BackendJobDescriptorKey(call, None, 2)
    val jobPathsAttempt = new JobPaths(wd, backendConfig, jobKeyAttempt, None)
    jobPathsAttempt.callRoot.toString shouldBe
      s"local-cromwell-executions/hello/$id/call-hello/attempt-2"

    val jobKeyShardedAttempt = new BackendJobDescriptorKey(call, Option(0), 2)
    val jobPathsShardedAttempt = new JobPaths(wd, backendConfig, jobKeyShardedAttempt, None)
    jobPathsShardedAttempt.callRoot.toString shouldBe
      s"local-cromwell-executions/hello/$id/call-hello/shard-0/attempt-2"
  }
}
