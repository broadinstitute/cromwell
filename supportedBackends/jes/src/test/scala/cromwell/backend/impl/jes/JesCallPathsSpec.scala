package cromwell.backend.impl.jes

import cromwell.backend.BackendSpec
import cromwell.util.SampleWdl
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.Mockito
import scala.concurrent.ExecutionContext.Implicits.global
import cromwell.backend.impl.jes.MockObjects._

class JesCallPathsSpec extends FlatSpec with Matchers with Mockito {

  import BackendSpec._
  import JesTestConfig._

  behavior of "JesCallPaths"

  it should "map the correct filenames" in {
    val workflowDescriptor = buildWorkflowDescriptor(SampleWdl.HelloWorld.wdlSource())
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)
    val jesConfiguration = new JesConfiguration(JesBackendConfigurationDescriptor)

    val callPaths = JesCallPaths(jobDescriptorKey, workflowDescriptor,
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

    val callPaths = JesCallPaths(jobDescriptorKey, workflowDescriptor, jesConfiguration)
    callPaths.returnCodePath.toUri.toString should
      be(s"gs://my-cromwell-workflows-bucket/hello/${workflowDescriptor.id}/call-hello/hello-rc.txt")
    callPaths.stdoutPath.toUri.toString should
      be(s"gs://my-cromwell-workflows-bucket/hello/${workflowDescriptor.id}/call-hello/hello-stdout.log")
    callPaths.stderrPath.toUri.toString should
      be(s"gs://my-cromwell-workflows-bucket/hello/${workflowDescriptor.id}/call-hello/hello-stderr.log")
    callPaths.jesLogPath.toUri.toString should
      be(s"gs://my-cromwell-workflows-bucket/hello/${workflowDescriptor.id}/call-hello/hello.log")
  }

  it should "map the correct call context" in {
    val workflowDescriptor = buildWorkflowDescriptor(SampleWdl.HelloWorld.wdlSource())
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)
    val jesConfiguration = new JesConfiguration(JesBackendConfigurationDescriptor)

    val callPaths = JesCallPaths(jobDescriptorKey, workflowDescriptor, jesConfiguration)
    callPaths.callContext.root.toUri.toString should
      be(s"gs://my-cromwell-workflows-bucket/hello/${workflowDescriptor.id}/call-hello")
    callPaths.callContext.stdout should be("hello-stdout.log")
    callPaths.callContext.stderr should be("hello-stderr.log")
  }

}
