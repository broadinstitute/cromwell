package cromwell.backend.google.pipelines.common

import com.google.cloud.NoCredentials
import cromwell.backend.BackendSpec
import cromwell.backend.io.JobPathsSpecHelper._
import cromwell.cloudsupport.gcp.auth.GoogleAuthModeSpec
import cromwell.core.TestKitSuite
import cromwell.util.SampleWdl
import org.scalatest.{FlatSpecLike, Matchers}
import org.specs2.mock.Mockito
import spray.json.{JsObject, JsString}

class PipelinesApiCallPathsSpec extends TestKitSuite with FlatSpecLike with Matchers with Mockito {

  import BackendSpec._
  import PipelinesApiTestConfig._

  behavior of "JesCallPaths"

  it should "map the correct filenames" in {
    GoogleAuthModeSpec.assumeHasApplicationDefaultCredentials()

    val workflowDescriptor = buildWdlWorkflowDescriptor(
      SampleWdl.HelloWorld.workflowSource(),
      inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.mapValues(JsString.apply)).compactPrint)
    )
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)
    val workflowPaths = PipelinesApiWorkflowPaths(workflowDescriptor, NoCredentials.getInstance(), NoCredentials.getInstance(), jesConfiguration, pathBuilders)
    
    val callPaths = PipelinesApiJobPaths(workflowPaths, jobDescriptorKey)
    
    callPaths.returnCodeFilename should be("rc")
    callPaths.stderr.getFileName.pathAsString should be("gs://my-cromwell-workflows-bucket/stderr")
    callPaths.stdout.getFileName.pathAsString should be("gs://my-cromwell-workflows-bucket/stdout")
    callPaths.jesLogFilename should be("hello.log")
  }

  it should "map the correct paths" in {
    GoogleAuthModeSpec.assumeHasApplicationDefaultCredentials()

    val workflowDescriptor = buildWdlWorkflowDescriptor(
      SampleWdl.HelloWorld.workflowSource(),
      inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.mapValues(JsString.apply)).compactPrint)
    )
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)
    val workflowPaths = PipelinesApiWorkflowPaths(workflowDescriptor, NoCredentials.getInstance(), NoCredentials.getInstance(), jesConfiguration, pathBuilders)

    val callPaths = PipelinesApiJobPaths(workflowPaths, jobDescriptorKey)
    
    callPaths.returnCode.pathAsString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/rc")
    callPaths.stdout.pathAsString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/stdout")
    callPaths.stderr.pathAsString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/stderr")
    callPaths.jesLogPath.pathAsString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/hello.log")
  }

  it should "map the correct call context" in {
    GoogleAuthModeSpec.assumeHasApplicationDefaultCredentials()

    val workflowDescriptor = buildWdlWorkflowDescriptor(
      SampleWdl.HelloWorld.workflowSource(),
      inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.mapValues(JsString.apply)).compactPrint)
    )
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)
    val workflowPaths = PipelinesApiWorkflowPaths(workflowDescriptor, NoCredentials.getInstance(), NoCredentials.getInstance(), jesConfiguration, pathBuilders)

    val callPaths = PipelinesApiJobPaths(workflowPaths, jobDescriptorKey)
    
    callPaths.callContext.root.pathAsString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello")
    callPaths.callContext.stdout should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/stdout")
    callPaths.callContext.stderr should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/stderr")
  }

}
