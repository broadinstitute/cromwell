package cromwell.backend.impl.jes

import cromwell.backend.BackendSpec
import cromwell.core.IntegrationTest
import cromwell.util.SampleWdl
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.Mockito


class JesCallPathsSpec extends FlatSpec with Matchers with Mockito with MockGcsFileSystemBuilder {

  import BackendSpec._
  import JesTestConfig._

  behavior of "JesCallPaths"

  it should "map the correct filenames" taggedAs IntegrationTest in {
    val workflowDescriptor = buildWorkflowDescriptor(SampleWdl.HelloWorld.wdlSource())
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)
    val jesConfiguration = new JesConfiguration(JesBackendConfigurationDescriptor)

    val callPaths = JesCallPaths(jobDescriptorKey, workflowDescriptor, jesConfiguration, buildMockedGcsFileSystem)
    callPaths.returnCodeFilename should be("hello-rc.txt")
    callPaths.stderrFilename should be("hello-stderr.log")
    callPaths.stdoutFilename should be("hello-stdout.log")
    callPaths.jesLogFilename should be("hello.log")
  }

  it should "map the correct paths" taggedAs IntegrationTest in {
    val workflowDescriptor = buildWorkflowDescriptor(SampleWdl.HelloWorld.wdlSource())
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)
    val jesConfiguration = new JesConfiguration(JesBackendConfigurationDescriptor)

    val callPaths = JesCallPaths(jobDescriptorKey, workflowDescriptor, jesConfiguration, buildMockedGcsFileSystem)
    callPaths.returnCodePath.toString should
      be(s"gs://my-cromwell-workflows-bucket/hello/${workflowDescriptor.id}/call-hello/hello-rc.txt")
    callPaths.stdoutPath.toString should
      be(s"gs://my-cromwell-workflows-bucket/hello/${workflowDescriptor.id}/call-hello/hello-stdout.log")
    callPaths.stderrPath.toString should
      be(s"gs://my-cromwell-workflows-bucket/hello/${workflowDescriptor.id}/call-hello/hello-stderr.log")
    callPaths.jesLogPath.toString should
      be(s"gs://my-cromwell-workflows-bucket/hello/${workflowDescriptor.id}/call-hello/hello.log")
  }

  it should "map the correct call context" taggedAs IntegrationTest in {
    val workflowDescriptor = buildWorkflowDescriptor(SampleWdl.HelloWorld.wdlSource())
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)
    val jesConfiguration = new JesConfiguration(JesBackendConfigurationDescriptor)

    val callPaths = JesCallPaths(jobDescriptorKey, workflowDescriptor, jesConfiguration, buildMockedGcsFileSystem)
    callPaths.callContext.root.toString should
      be(s"gs://my-cromwell-workflows-bucket/hello/${workflowDescriptor.id}/call-hello")
    callPaths.callContext.stdout should be("hello-stdout.log")
    callPaths.callContext.stderr should be("hello-stderr.log")
  }

}
