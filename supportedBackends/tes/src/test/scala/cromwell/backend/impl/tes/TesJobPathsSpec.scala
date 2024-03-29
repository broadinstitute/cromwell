package cromwell.backend.impl.tes

import better.files._
import common.assertion.CromwellTimeoutSpec
import cromwell.backend.io.JobPathsSpecHelper._
import cromwell.backend.{BackendJobDescriptorKey, BackendSpec}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wom.graph.CommandCallNode

class TesJobPathsSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with BackendSpec {

  "JobPaths" should "provide correct paths for a job" in {

    val wd = buildWdlWorkflowDescriptor(TestWorkflows.HelloWorld)
    val call: CommandCallNode = wd.callable.taskCallNodes.head
    val jobKey = BackendJobDescriptorKey(call, None, 1)
    val jobPaths = TesJobPaths(jobKey, wd, TesTestConfig.backendConfig)
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
    jobPaths.azureLogPathsForMetadata shouldBe
      Map(
        "tes_stdout" -> File(s"local-cromwell-executions/wf_hello/$id/call-hello/tes_task/stdout.txt").pathAsString,
        "tes_stderr" -> File(s"local-cromwell-executions/wf_hello/$id/call-hello/tes_task/stderr.txt").pathAsString
      )

    val jobKeySharded = BackendJobDescriptorKey(call, Option(0), 1)
    val jobPathsSharded = TesJobPaths(jobKeySharded, wd, TesTestConfig.backendConfig)
    jobPathsSharded.callExecutionRoot.toString shouldBe
      File(s"local-cromwell-executions/wf_hello/$id/call-hello/shard-0/execution").pathAsString

    val jobKeyAttempt = BackendJobDescriptorKey(call, None, 2)
    val jobPathsAttempt = TesJobPaths(jobKeyAttempt, wd, TesTestConfig.backendConfig)
    jobPathsAttempt.callExecutionRoot.toString shouldBe
      File(s"local-cromwell-executions/wf_hello/$id/call-hello/attempt-2/execution").pathAsString

    val jobKeyShardedAttempt = BackendJobDescriptorKey(call, Option(0), 2)
    val jobPathsShardedAttempt = TesJobPaths(jobKeyShardedAttempt, wd, TesTestConfig.backendConfig)
    jobPathsShardedAttempt.callExecutionRoot.toString shouldBe
      File(s"local-cromwell-executions/wf_hello/$id/call-hello/shard-0/attempt-2/execution").pathAsString
  }
}
