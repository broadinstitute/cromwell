package cromwell.backend.google.pipelines.batch.models

import com.google.cloud.NoCredentials
import common.collections.EnhancedCollections._
import cromwell.backend.BackendSpec
import cromwell.backend.google.pipelines.batch.actors.GcpBatchInitializationActor
import cromwell.backend.google.pipelines.batch.models.GcpBatchTestConfig.{batchConfiguration, pathBuilders}
import cromwell.backend.io.JobPathsSpecHelper._
import cromwell.core.TestKitSuite
import cromwell.util.SampleWdl
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import spray.json.{JsObject, JsString}

import scala.concurrent.ExecutionContext.Implicits.global

class GcpBatchJobPathsSpec extends TestKitSuite with AnyFlatSpecLike with Matchers {

  import BackendSpec._

  behavior of "JesCallPaths"

  it should "map the correct filenames" in (pending) {
    val workflowDescriptor = buildWdlWorkflowDescriptor(
      SampleWdl.HelloWorld.workflowSource(),
      inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.safeMapValues(JsString.apply)).compactPrint)
    )
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)

    // TODO: review if the method defaultStandardStreamNameToFileNameMetadataMapper is required in the initialization actor
    val workflowPaths = GcpBatchWorkflowPaths(workflowDescriptor, NoCredentials.getInstance(), NoCredentials.getInstance(), batchConfiguration, pathBuilders(), GcpBatchInitializationActor.defaultStandardStreamNameToFileNameMetadataMapper)

    val callPaths = GcpBatchJobPaths(workflowPaths, jobDescriptorKey)

    callPaths.returnCodeFilename should be("rc")
    callPaths.stderr.getFileName.pathAsString should be("gs://my-cromwell-workflows-bucket/stderr")
    callPaths.stdout.getFileName.pathAsString should be("gs://my-cromwell-workflows-bucket/stdout")
    callPaths.jesLogFilename should be("hello.log")
  }

  it should "map the correct paths" in (pending) {
    val workflowDescriptor = buildWdlWorkflowDescriptor(
      SampleWdl.HelloWorld.workflowSource(),
      inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.safeMapValues(JsString.apply)).compactPrint)
    )
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)

    // TODO: review if the method defaultStandardStreamNameToFileNameMetadataMapper is required in the initialization actor
    val workflowPaths = GcpBatchWorkflowPaths(workflowDescriptor, NoCredentials.getInstance(), NoCredentials.getInstance(), batchConfiguration, pathBuilders(), GcpBatchInitializationActor.defaultStandardStreamNameToFileNameMetadataMapper)

    val callPaths = GcpBatchJobPaths(workflowPaths, jobDescriptorKey)

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
    pending
    val workflowDescriptor = buildWdlWorkflowDescriptor(
      SampleWdl.HelloWorld.workflowSource(),
      inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.safeMapValues(JsString.apply)).compactPrint)
    )
    val jobDescriptorKey = firstJobDescriptorKey(workflowDescriptor)

    // TODO: review if the method defaultStandardStreamNameToFileNameMetadataMapper is required in the initialization actor
    val workflowPaths = GcpBatchWorkflowPaths(workflowDescriptor, NoCredentials.getInstance(), NoCredentials.getInstance(), batchConfiguration, pathBuilders(), GcpBatchInitializationActor.defaultStandardStreamNameToFileNameMetadataMapper)

    val callPaths = GcpBatchJobPaths(workflowPaths, jobDescriptorKey)

    callPaths.callContext.root.pathAsString should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello")
    callPaths.callContext.stdout should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/stdout")
    callPaths.callContext.stderr should
      be(s"gs://my-cromwell-workflows-bucket/wf_hello/${workflowDescriptor.id}/call-hello/stderr")
  }

}
