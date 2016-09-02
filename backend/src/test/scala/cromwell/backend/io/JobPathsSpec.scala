package cromwell.backend.io

import java.nio.file.Paths

import better.files._
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
  val defaultBackendConfigDescriptor = BackendConfigurationDescriptor(backendConfig, globalConfig)

  "JobPaths" should "provide correct paths for a job" in {

    val wd = buildWorkflowDescriptor(TestWorkflows.HelloWorld)
    val call: Call = wd.workflowNamespace.workflow.calls.head
    val jobKey = BackendJobDescriptorKey(call, None, 1)
    val jobPaths = new JobPaths(wd, backendConfig, jobKey)
    val id = wd.id
    jobPaths.callRoot.toString shouldBe
      File(s"local-cromwell-executions/hello/$id/call-hello").fullPath
    jobPaths.returnCode.toString shouldBe
      File(s"local-cromwell-executions/hello/$id/call-hello/rc").fullPath
    jobPaths.script.toString shouldBe
      File(s"local-cromwell-executions/hello/$id/call-hello/script").fullPath
    jobPaths.stderr.toString shouldBe
      File(s"local-cromwell-executions/hello/$id/call-hello/stderr").fullPath
    jobPaths.stdout.toString shouldBe
      File(s"local-cromwell-executions/hello/$id/call-hello/stdout").fullPath
    jobPaths.callRoot.toString shouldBe
      File(s"local-cromwell-executions/hello/$id/call-hello").fullPath
    jobPaths.callDockerRoot.toString shouldBe
      File(s"/root/hello/$id/call-hello").fullPath
    jobPaths.toDockerPath(Paths.get(s"local-cromwell-executions/hello/$id/call-hello/stdout")).toString shouldBe
      File(s"/root/hello/$id/call-hello/stdout").fullPath
    jobPaths.toDockerPath(Paths.get("/root/dock/path")).toString shouldBe
      File("/root/dock/path").fullPath

    val jobKeySharded = BackendJobDescriptorKey(call, Option(0), 1)
    val jobPathsSharded = new JobPaths(wd, backendConfig, jobKeySharded)
    jobPathsSharded.callRoot.toString shouldBe
      File(s"local-cromwell-executions/hello/$id/call-hello/shard-0").fullPath

    val jobKeyAttempt = BackendJobDescriptorKey(call, None, 2)
    val jobPathsAttempt = new JobPaths(wd, backendConfig, jobKeyAttempt)
    jobPathsAttempt.callRoot.toString shouldBe
      File(s"local-cromwell-executions/hello/$id/call-hello/attempt-2").fullPath

    val jobKeyShardedAttempt = BackendJobDescriptorKey(call, Option(0), 2)
    val jobPathsShardedAttempt = new JobPaths(wd, backendConfig, jobKeyShardedAttempt)
    jobPathsShardedAttempt.callRoot.toString shouldBe
      File(s"local-cromwell-executions/hello/$id/call-hello/shard-0/attempt-2").fullPath
  }
}
