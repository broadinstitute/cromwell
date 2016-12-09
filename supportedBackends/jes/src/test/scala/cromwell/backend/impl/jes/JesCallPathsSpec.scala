package cromwell.backend.impl.jes

import cromwell.backend.BackendSpec
import cromwell.core.TestKitSuite
import cromwell.core.path.PathImplicits._
import cromwell.util.SampleWdl
import org.scalatest.{FlatSpecLike, Matchers}
import org.specs2.mock.Mockito

class JesCallPathsSpec extends TestKitSuite with FlatSpecLike with Matchers with Mockito {

  import BackendSpec._
  import JesTestConfig._

  behavior of "JesCallPaths"

  it should "map the correct filenames" in {
    val workflowDescriptor = buildWorkflowDescriptor(SampleWdl.HelloWorld.wdlSource())
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)
    val jesConfiguration = new JesConfiguration(JesBackendConfigurationDescriptor)

    val callPaths = JesJobPaths(jobDescriptorKey, workflowDescriptor,
      jesConfiguration)
    callPaths.returnCodeFilename should be("hello-rc.txt")
    callPaths.stderrFilename should be("hello-stderr.log")
    callPaths.stdoutFilename should be("hello-stdout.log")
    callPaths.jesLogFilename should be("hello.log")
  }

  it should "map the correct paths" in {
    val workflowDescriptor = buildWorkflowDescriptor(SampleWdl.HelloWorld.wdlSource())
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)
    val jesConfiguration = new JesConfiguration(JesBackendConfigurationDescriptor)

    val callPaths = JesJobPaths(jobDescriptorKey, workflowDescriptor, jesConfiguration)
    callPaths.returnCode.toRealString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/hello-rc.txt")
    callPaths.stdout.toRealString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/hello-stdout.log")
    callPaths.stderr.toRealString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/hello-stderr.log")
    callPaths.jesLogPath.toRealString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/hello.log")
  }

  it should "map the correct call context" in {
    val workflowDescriptor = buildWorkflowDescriptor(SampleWdl.HelloWorld.wdlSource())
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)
    val jesConfiguration = new JesConfiguration(JesBackendConfigurationDescriptor)

    val callPaths = JesJobPaths(jobDescriptorKey, workflowDescriptor, jesConfiguration)
    callPaths.callContext.root.toRealString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello")
    callPaths.callContext.stdout should be("hello-stdout.log")
    callPaths.callContext.stderr should be("hello-stderr.log")
  }

}
