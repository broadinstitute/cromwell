package cromwell.backend.impl.jes

import com.google.cloud.NoCredentials
import cromwell.backend.BackendSpec
import cromwell.cloudsupport.gcp.auth.GoogleAuthModeSpec
import cromwell.core.TestKitSuite
import cromwell.util.SampleWdl
import org.scalatest.{FlatSpecLike, Matchers}
import org.specs2.mock.Mockito

class JesCallPathsSpec extends TestKitSuite with FlatSpecLike with Matchers with Mockito {

  import BackendSpec._
  import JesTestConfig._

  behavior of "JesCallPaths"

  it should "map the correct filenames" in {
    GoogleAuthModeSpec.assumeHasApplicationDefaultCredentials()

    val workflowDescriptor = buildWorkflowDescriptor(SampleWdl.HelloWorld.workflowSource())
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)
    val jesConfiguration = new JesConfiguration(JesBackendConfigurationDescriptor)
    val workflowPaths = JesWorkflowPaths(workflowDescriptor, NoCredentials.getInstance(), NoCredentials.getInstance(), jesConfiguration)
    
    val callPaths = JesJobPaths(workflowPaths, jobDescriptorKey)
    
    callPaths.returnCodeFilename should be("hello-rc.txt")
    callPaths.stderrFilename should be("hello-stderr.log")
    callPaths.stdoutFilename should be("hello-stdout.log")
    callPaths.jesLogFilename should be("hello.log")
  }

  it should "map the correct paths" in {
    GoogleAuthModeSpec.assumeHasApplicationDefaultCredentials()

    val workflowDescriptor = buildWorkflowDescriptor(SampleWdl.HelloWorld.workflowSource())
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)
    val jesConfiguration = new JesConfiguration(JesBackendConfigurationDescriptor)
    val workflowPaths = JesWorkflowPaths(workflowDescriptor, NoCredentials.getInstance(), NoCredentials.getInstance(), jesConfiguration)

    val callPaths = JesJobPaths(workflowPaths, jobDescriptorKey)
    
    callPaths.returnCode.pathAsString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/hello-rc.txt")
    callPaths.stdout.pathAsString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/hello-stdout.log")
    callPaths.stderr.pathAsString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/hello-stderr.log")
    callPaths.jesLogPath.pathAsString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/hello.log")
  }

  it should "map the correct call context" in {
    GoogleAuthModeSpec.assumeHasApplicationDefaultCredentials()

    val workflowDescriptor = buildWorkflowDescriptor(SampleWdl.HelloWorld.workflowSource())
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)
    val jesConfiguration = new JesConfiguration(JesBackendConfigurationDescriptor)
    val workflowPaths = JesWorkflowPaths(workflowDescriptor, NoCredentials.getInstance(), NoCredentials.getInstance(), jesConfiguration)

    val callPaths = JesJobPaths(workflowPaths, jobDescriptorKey)
    
    callPaths.callContext.root.pathAsString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello")
    callPaths.callContext.stdout should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/hello-stdout.log")
    callPaths.callContext.stderr should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/hello-stderr.log")
  }

}
