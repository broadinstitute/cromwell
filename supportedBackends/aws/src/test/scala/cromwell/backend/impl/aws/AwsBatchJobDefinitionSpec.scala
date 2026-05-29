package cromwell.backend.impl.aws

import common.mock.MockSugar
import cromwell.backend.BackendSpec._
import cromwell.backend.impl.aws.io.AwsBatchWorkingDisk
import cromwell.backend.io.JobPaths
import cromwell.backend.{BackendJobDescriptor, BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core.callcaching.NoDocker
import cromwell.util.SampleWdl
import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.refineMV
import org.scalatest.matchers.should.Matchers
import org.scalatest.wordspec.AnyWordSpecLike
import spray.json.{JsObject, JsString}
import wdl4s.parser.MemoryUnit
import wom.format.MemorySize
import wom.graph.CommandCallNode

import scala.jdk.CollectionConverters._

class AwsBatchJobDefinitionSpec extends AnyWordSpecLike with Matchers with MockSugar {

  // Create a simple workflow descriptor for testing
  val workflowDescriptor: BackendWorkflowDescriptor = buildWdlWorkflowDescriptor(
    SampleWdl.HelloWorld.workflowSource(),
    inputFileAsJson = Option(JsObject(SampleWdl.HelloWorld.rawInputs.map { case (k, v) =>
      k -> JsString(v)
    }).compactPrint)
  )

  val call: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
  val jobKey: BackendJobDescriptorKey = BackendJobDescriptorKey(call, None, 1)

  // Mock job descriptor for testing
  val jobDescriptor: BackendJobDescriptor = BackendJobDescriptor(
    workflowDescriptor,
    jobKey,
    Map.empty,
    Map.empty,
    NoDocker,
    None,
    Map.empty
  )

  // Test runtime attributes
  val runtimeAttributes = AwsBatchRuntimeAttributes(
    cpu = refineMV[Positive](1),
    gpuCount = 0,
    zones = Vector("us-east-1a"),
    memory = MemorySize(2.0, MemoryUnit.GB),
    disks = Seq.empty,
    dockerImage = "ubuntu:latest",
    queueArn = "arn:aws:batch:us-east-1:123456789:job-queue/test-queue",
    failOnStderr = false,
    continueOnReturnCode = cromwell.backend.validation.ContinueOnReturnCodeSet(Set(0)),
    noAddress = false,
    scriptS3BucketName = "test-bucket",
    awsBatchRetryAttempts = 1,
    awsBatchEvaluateOnExit = Vector(Map.empty[String, String]),
    ulimits = Vector(Map.empty[String, String]),
    efsDelocalize = false,
    efsMakeMD5 = false,
    fuseMount = false,
    sharedMemorySize = MemorySize(10, MemoryUnit.MB),
    jobTimeout = 0,
    logGroupName = "/aws/batch/job",
    additionalTags = Map.empty,
    fileSystem = "s3"
  )

  // Mock job paths for testing
  val mockJobPaths: JobPaths = mock[JobPaths]

  private def buildContext(attrs: AwsBatchRuntimeAttributes): AwsBatchJobDefinitionContext =
    AwsBatchJobDefinitionContext(
      runtimeAttributes = attrs,
      commandText = "echo hello",
      dockerRcPath = "/tmp/rc.txt",
      dockerStdoutPath = "/tmp/stdout.log",
      dockerStderrPath = "/tmp/stderr.log",
      jobDescriptor = jobDescriptor,
      jobPaths = mockJobPaths,
      inputs = Set.empty,
      outputs = Set.empty,
      fsxMntPoint = None,
      None,
      None,
      None,
      None,
      None,
      None
    )

  "AwsBatchJobDefinition StandardAwsBatchJobDefinitionBuilder" should {

    "build a job definition name without jobRoleArn when not provided" in {

      val context = AwsBatchJobDefinitionContext(
        runtimeAttributes = runtimeAttributes,
        commandText = "echo hello",
        dockerRcPath = "/tmp/rc.txt",
        dockerStdoutPath = "/tmp/stdout.log",
        dockerStderrPath = "/tmp/stderr.log",
        jobDescriptor = jobDescriptor,
        jobPaths = mockJobPaths,
        inputs = Set.empty,
        outputs = Set.empty,
        fsxMntPoint = None,
        None,
        None,
        None,
        None,
        None,
        None
      )

      val jobDefinition = StandardAwsBatchJobDefinitionBuilder.build(context)

      // Job definition name should be deterministic based on inputs
      // The colon in "ubuntu:latest" gets replaced with underscore
      jobDefinition.name should startWith("cromwell_ubuntu_latest")
      // The name should contain the SHA1 hash at the end
      jobDefinition.name should fullyMatch regex """cromwell_ubuntu_latest_[a-f0-9]{40}"""
    }

    "build a job definition name with jobRoleArn when provided" in {
      // This test makes sure that the Job Definition is rebuilt correctly
      // if the user submits a workflow with different JobRoleARN.
      val roleArn = "arn:aws:iam::123456789012:role/MyJobRole"

      val context = AwsBatchJobDefinitionContext(
        runtimeAttributes = runtimeAttributes,
        commandText = "echo hello",
        dockerRcPath = "/tmp/rc.txt",
        dockerStdoutPath = "/tmp/stdout.log",
        dockerStderrPath = "/tmp/stderr.log",
        jobDescriptor = jobDescriptor,
        jobPaths = mockJobPaths,
        inputs = Set.empty,
        outputs = Set.empty,
        fsxMntPoint = None,
        None,
        None,
        None,
        None,
        None,
        Some(roleArn)
      )

      val jobDefinitionWithRole = StandardAwsBatchJobDefinitionBuilder.build(context)

      // Job definition name should be different when roleArn is included
      // The colon in "ubuntu:latest" gets replaced with underscore
      jobDefinitionWithRole.name should startWith("cromwell_ubuntu_latest")
      // The name should contain the SHA1 hash at the end
      jobDefinitionWithRole.name should fullyMatch regex """cromwell_ubuntu_latest_[a-f0-9]{40}"""

      // Create the same job definition without roleArn to compare
      val contextNoRole = AwsBatchJobDefinitionContext(
        runtimeAttributes = runtimeAttributes,
        commandText = "echo hello",
        dockerRcPath = "/tmp/rc.txt",
        dockerStdoutPath = "/tmp/stdout.log",
        dockerStderrPath = "/tmp/stderr.log",
        jobDescriptor = jobDescriptor,
        jobPaths = mockJobPaths,
        inputs = Set.empty,
        outputs = Set.empty,
        fsxMntPoint = None,
        None,
        None,
        None,
        None,
        None,
        None
      )

      val jobDefinitionNoRole = StandardAwsBatchJobDefinitionBuilder.build(contextNoRole)

      // The two names should be different because roleArn is included in the hash
      jobDefinitionWithRole.name should not equal jobDefinitionNoRole.name
    }

    "apply jobRoleArn to container properties when provided" in {
      val roleArn = "arn:aws:iam::123456789012:role/MyJobRole"

      val context = AwsBatchJobDefinitionContext(
        runtimeAttributes = runtimeAttributes,
        commandText = "echo hello",
        dockerRcPath = "/tmp/rc.txt",
        dockerStdoutPath = "/tmp/stdout.log",
        dockerStderrPath = "/tmp/stderr.log",
        jobDescriptor = jobDescriptor,
        jobPaths = mockJobPaths,
        inputs = Set.empty,
        outputs = Set.empty,
        fsxMntPoint = None,
        None,
        None,
        None,
        None,
        None,
        Some(roleArn)
      )

      val jobDefinition = StandardAwsBatchJobDefinitionBuilder.build(context)
      val containerProperties = jobDefinition.containerProperties

      // Verify the job role ARN is set on the container properties
      containerProperties.jobRoleArn() should equal(roleArn)
    }

    "not set jobRoleArn on container properties when not provided" in {
      val context = AwsBatchJobDefinitionContext(
        runtimeAttributes = runtimeAttributes,
        commandText = "echo hello",
        dockerRcPath = "/tmp/rc.txt",
        dockerStdoutPath = "/tmp/stdout.log",
        dockerStderrPath = "/tmp/stderr.log",
        jobDescriptor = jobDescriptor,
        jobPaths = mockJobPaths,
        inputs = Set.empty,
        outputs = Set.empty,
        fsxMntPoint = None,
        None,
        None,
        None,
        None,
        None,
        None
      )

      val jobDefinition = StandardAwsBatchJobDefinitionBuilder.build(context)
      val containerProperties = jobDefinition.containerProperties

      // Verify no job role ARN is set when not provided
      Option(containerProperties.jobRoleArn()).isEmpty should be(true)
    }

    "set CROMWELL_DISK_GB environment variable when working disk sizeGb > 0" in {
      val context = buildContext(runtimeAttributes.copy(disks = Seq(AwsBatchWorkingDisk(sizeGb = Some(500)))))
      val containerProperties = StandardAwsBatchJobDefinitionBuilder.build(context).containerProperties
      val envNames = containerProperties.environment().asScala.map(_.name())

      envNames should contain("CROMWELL_DISK_GB")
      containerProperties
        .environment()
        .asScala
        .find(_.name() == "CROMWELL_DISK_GB")
        .map(_.value()) shouldEqual Some("500")
    }

    "not set CROMWELL_DISK_GB environment variable when working disk sizeGb is 0" in {
      val context = buildContext(runtimeAttributes.copy(disks = Seq(AwsBatchWorkingDisk())))
      val containerProperties = StandardAwsBatchJobDefinitionBuilder.build(context).containerProperties

      containerProperties.environment().asScala.map(_.name()) should not contain "CROMWELL_DISK_GB"
    }

    "mount cromwellDiskUtils volume and set privileged when working disk sizeGb > 0" in {
      val context = buildContext(runtimeAttributes.copy(disks = Seq(AwsBatchWorkingDisk(sizeGb = Some(500)))))
      val containerProperties = StandardAwsBatchJobDefinitionBuilder.build(context).containerProperties

      containerProperties.volumes().asScala.map(_.name()) should contain("cromwellDiskUtils")
      containerProperties.mountPoints().asScala.map(_.sourceVolume()) should contain("cromwellDiskUtils")
      containerProperties.mountPoints().asScala.map(_.containerPath()) should contain("/usr/local/cromwell-disk-utils")
      containerProperties.privileged() shouldBe true
    }

    "not mount cromwellDiskUtils volume and not set privileged when working disk sizeGb is 0" in {
      val context = buildContext(runtimeAttributes)
      val containerProperties = StandardAwsBatchJobDefinitionBuilder.build(context).containerProperties

      containerProperties.volumes().asScala.map(_.name()) should not contain "cromwellDiskUtils"
      containerProperties.mountPoints().asScala.map(_.sourceVolume()) should not contain "cromwellDiskUtils"
      containerProperties.privileged() shouldBe false
    }

    "produce distinct job definition names for different working disk sizes" in {
      val context100 = buildContext(runtimeAttributes.copy(disks = Seq(AwsBatchWorkingDisk(sizeGb = Some(100)))))
      val context500 = buildContext(runtimeAttributes.copy(disks = Seq(AwsBatchWorkingDisk(sizeGb = Some(500)))))

      val name100 = StandardAwsBatchJobDefinitionBuilder.build(context100).name
      val name500 = StandardAwsBatchJobDefinitionBuilder.build(context500).name

      name100 should not equal name500
    }
  }
}
