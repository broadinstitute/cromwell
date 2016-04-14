package cromwell.backend.impl.local

import java.nio.file.Paths

import cromwell.backend.BackendJobDescriptorKey
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.Call

class JobPathsSpec extends FlatSpec with Matchers with BackendTestkitSpec {

  "JobPaths" should "provide correct paths for a job" in {

    val wd = buildWorkflowDescriptor(TestWorkflows.HelloWorld)
    val call: Call = wd.workflowNamespace.workflow.calls.head
    val jobKey = new BackendJobDescriptorKey(call, None, 1)
    val jobPaths = new JobPaths(wd, defaultConfig, jobKey)
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
    val jobPathsSharded = new JobPaths(wd, defaultConfig, jobKeySharded)
    jobPathsSharded.callRoot.toString shouldBe
      s"local-cromwell-executions/hello/$id/call-hello/shard-0"

    val jobKeyAttempt = new BackendJobDescriptorKey(call, None, 2)
    val jobPathsAttempt = new JobPaths(wd, defaultConfig, jobKeyAttempt)
    jobPathsAttempt.callRoot.toString shouldBe
      s"local-cromwell-executions/hello/$id/call-hello/attempt-2"

    val jobKeyShardedAttempt = new BackendJobDescriptorKey(call, Option(0), 2)
    val jobPathsShardedAttempt = new JobPaths(wd, defaultConfig, jobKeyShardedAttempt)
    jobPathsShardedAttempt.callRoot.toString shouldBe
      s"local-cromwell-executions/hello/$id/call-hello/shard-0/attempt-2"
  }
}
