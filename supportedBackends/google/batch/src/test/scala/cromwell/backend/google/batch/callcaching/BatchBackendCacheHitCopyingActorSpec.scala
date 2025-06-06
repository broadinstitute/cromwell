package cromwell.backend.google.batch.callcaching

import akka.event.NoLogging
import akka.testkit.{ImplicitSender, TestFSMRef, TestProbe}
import com.typesafe.config.ConfigFactory
import common.mock.MockSugar
import cromwell.backend.BackendCacheHitCopyingActor.{CopyingOutputsFailedResponse, CopyOutputsCommand}
import cromwell.backend.BackendJobExecutionActor.JobSucceededResponse
import cromwell.backend.google.batch.models.GcpBatchConfigurationAttributes.VirtualPrivateCloudConfiguration
import cromwell.backend.google.batch.models._
import cromwell.backend.io.JobPaths
import cromwell.backend.standard.StandardValidatedRuntimeAttributesBuilder
import cromwell.backend.standard.callcaching.CopyingActorBlacklistCacheSupport.HasFormatting
import cromwell.backend.standard.callcaching.StandardCacheHitCopyingActor._
import cromwell.backend.standard.callcaching._
import cromwell.backend.validation.ValidatedRuntimeAttributes
import cromwell.backend.{BackendJobDescriptor, BackendJobDescriptorKey, BackendWorkflowDescriptor}
import cromwell.core._
import cromwell.core.callcaching.DockerWithHash
import cromwell.core.io.DefaultIoCommand.DefaultIoCopyCommand
import cromwell.core.io.{IoFailure, IoReadForbiddenFailure, IoSuccess}
import cromwell.core.path.Path
import cromwell.services.CallCaching.CallCachingEntryId
import cromwell.services.instrumentation.CromwellCount
import cromwell.services.instrumentation.InstrumentationService.InstrumentationServiceMessage
import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.refineMV
import org.mockito.ArgumentMatchers._
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.slf4j.Logger
import wom.callable.CommandTaskDefinition
import wom.graph.{CommandCallNode, FullyQualifiedName, LocalName, WomIdentifier}
import wom.values.WomValue

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Success, Try}

class BatchBackendCacheHitCopyingActorSpec
    extends TestKitSuite
    with AnyFlatSpecLike
    with Matchers
    with ImplicitSender
    with MockSugar
    with Eventually {

  behavior of "BatchBackendCacheHitCopyingActor"

  private val LockedDownBucket = "locked-down-bucket"
  private val WideOpenBucket = "wide-open-bucket"
  private val GoogleProject = "cache_as_cache_can"

  it should "do all the right things with blacklisting hits and buckets with groupings enabled" in {

    val configString =
      """
        |call-caching {
        |  enabled: true
        |  blacklist-cache {
        |    enabled: true
        |
        |    groupings {
        |      workflow-option: google_project
        |    }
        |  }
        |}
        |""".stripMargin

    val blacklistManager = new CallCachingBlacklistManager(ConfigFactory.parseString(configString), NoLogging)
    val grouping = Option(GoogleProject)
    val workflow = buildWorkflow(grouping)
    val blacklistCache = blacklistManager.blacklistCacheFor(workflow).get

    // Make sure we got the expected type of cache
    blacklistCache match {
      case _: GroupingBlacklistCache =>
      case bad =>
        fail(s"Unexpected blacklist cache type, expected GroupingBlacklistCache: ${bad.getClass.getSimpleName}")
    }

    {
      // Step 0: a successful copy attempt. There's a lot of darkness ahead so begin with a bit of light.
      val ioActor = TestProbe("ioActor")
      val serviceRegistryActor = TestProbe("serviceRegistryActor")
      val supervisor = TestProbe("supervisor")
      val copyActor = buildCopyActor(
        workflow = workflow,
        blacklistCache = blacklistCache,
        fakeIoActor = ioActor,
        fakeServiceRegistryActor = serviceRegistryActor,
        supervisor = supervisor,
        grouping = grouping
      )

      val copyCommand = buildCopyCommand(hitId = 0, bucket = WideOpenBucket)
      supervisor watch copyActor

      copyActor ! copyCommand

      eventually {
        copyActor.underlyingActor.stateName shouldBe WaitingForIoResponses
      }

      ioActor.expectMsgPF(5 seconds) { case ioCommand: DefaultIoCopyCommand =>
        ioActor.reply(IoSuccess(ioCommand, ()))
      }

      supervisor.expectMsgPF(5 seconds) { case _: JobSucceededResponse => }

      val counts = instrumentationCounts(n = 4, serviceRegistryActor = serviceRegistryActor)
      val (List(hitBegin, hitEnd), List(bucketBegin, bucketEnd)) = counts partition {
        _.bucket.path.toList.contains("hit")
      }

      hitBegin.bucket.path.toList shouldBe expectedMetric(Hit, Read, UntestedCacheResult)
      bucketBegin.bucket.path.toList shouldBe expectedMetric(Bucket, Read, UntestedCacheResult)

      hitEnd.bucket.path.toList shouldBe expectedMetric(Hit, Write, GoodCacheResult)
      bucketEnd.bucket.path.toList shouldBe expectedMetric(Bucket, Write, GoodCacheResult)

      blacklistCache.bucketCache.size() shouldBe 1
      blacklistCache.bucketCache.get(WideOpenBucket) shouldBe GoodCacheResult

      blacklistCache.hitCache.size() shouldBe 1
      blacklistCache.hitCache.get(CallCachingEntryId(0)) shouldBe GoodCacheResult
    }

    {
      // Step 1: an attempt to read a hit of unknown status from a bucket of unknown status, but the IoActor will report
      // a forbidden (403) failure which should cause hit and bucket blacklisting.
      val ioActor = TestProbe("ioActor")
      val serviceRegistryActor = TestProbe("serviceRegistryActor")
      val supervisor = TestProbe("supervisor")

      val copyActor = buildCopyActor(
        workflow = workflow,
        blacklistCache = blacklistCache,
        fakeIoActor = ioActor,
        fakeServiceRegistryActor = serviceRegistryActor,
        supervisor = supervisor,
        grouping = grouping
      )

      val command = buildCopyCommand(hitId = 1, bucket = LockedDownBucket)
      supervisor watch copyActor

      copyActor ! command

      eventually {
        copyActor.underlyingActor.stateName shouldBe WaitingForIoResponses
      }

      ioActor.expectMsgPF(5 seconds) { case ioCommand: DefaultIoCopyCommand =>
        val failedPath = command.jobDetritusFiles(JobPaths.ReturnCodePathKey)
        ioActor.reply(IoReadForbiddenFailure(ioCommand, new RuntimeException(), failedPath))
      }

      supervisor.expectMsgPF(5 seconds) { case _: CopyingOutputsFailedResponse => }

      val counts = instrumentationCounts(n = 4, serviceRegistryActor = serviceRegistryActor)

      // Expect read hit and read bucket UntestedCacheResult followed by write hit and write bucket BadCacheResult.
      {
        val (List(hitBegin, hitEnd), List(bucketBegin, bucketEnd)) = counts partition {
          _.bucket.path.toList.contains("hit")
        }

        hitBegin.bucket.path.toList shouldBe expectedMetric(Hit, Read, UntestedCacheResult)
        bucketBegin.bucket.path.toList shouldBe expectedMetric(Bucket, Read, UntestedCacheResult)

        hitEnd.bucket.path.toList shouldBe expectedMetric(Hit, Write, BadCacheResult)
        bucketEnd.bucket.path.toList shouldBe expectedMetric(Bucket, Write, BadCacheResult)
      }

      // Assert blacklist entries were made for bucket and hit.
      blacklistCache.bucketCache.size() shouldBe 2
      blacklistCache.bucketCache.get(WideOpenBucket) shouldBe GoodCacheResult
      blacklistCache.bucketCache.get(LockedDownBucket) shouldBe BadCacheResult

      blacklistCache.hitCache.size() shouldBe 2
      blacklistCache.hitCache.get(CallCachingEntryId(0)) shouldBe GoodCacheResult
      blacklistCache.hitCache.get(CallCachingEntryId(1)) shouldBe BadCacheResult

      supervisor.expectTerminated(copyActor, 5 seconds)
    }

    {
      // Step 2: an attempt to read an unknown hit from a blacklisted bucket.
      val ioActor = TestProbe("ioActor")
      val serviceRegistryActor = TestProbe("serviceRegistryActor")
      val supervisor = TestProbe("supervisor")

      val copyActor = buildCopyActor(
        workflow = workflow,
        blacklistCache = blacklistCache,
        fakeIoActor = ioActor,
        fakeServiceRegistryActor = serviceRegistryActor,
        supervisor = supervisor,
        grouping = grouping
      )

      supervisor watch copyActor

      val command = buildCopyCommand(hitId = 2, bucket = LockedDownBucket)
      copyActor ! command

      supervisor.expectMsgPF(5 seconds) { case _: CopyingOutputsFailedResponse => }
      // In this circumstance the copy actor just stops itself without transitioning out of Idle.
      supervisor.expectTerminated(copyActor)
      // Copying should be short-circuited by the bucket being blacklisted, so no communication with the IoActor.
      ioActor.expectNoMessage(max = 5 seconds)

      val List(hitMessage, bucketMessage) = instrumentationCounts(n = 2, serviceRegistryActor = serviceRegistryActor)

      // Hit status is unknown but bucket status is known bad.
      hitMessage.bucket.path.toList shouldBe expectedMetric(Hit, Read, UntestedCacheResult)
      bucketMessage.bucket.path.toList shouldBe expectedMetric(Bucket, Read, BadCacheResult)

      blacklistCache.bucketCache.size() shouldBe 2
      blacklistCache.bucketCache.get(WideOpenBucket) shouldBe GoodCacheResult
      blacklistCache.bucketCache.get(LockedDownBucket) shouldBe BadCacheResult

      blacklistCache.hitCache.size() shouldBe 3
      blacklistCache.hitCache.get(CallCachingEntryId(0)) shouldBe GoodCacheResult
      blacklistCache.hitCache.get(CallCachingEntryId(1)) shouldBe BadCacheResult
      blacklistCache.hitCache.get(CallCachingEntryId(2)) shouldBe UntestedCacheResult
    }

    {
      // Step 3: a generic failure to read a cache hit from a bucket not known to be bad should cause the hit to be
      // marked bad but not its containing bucket.
      val ioActor = TestProbe("ioActor")
      val serviceRegistryActor = TestProbe("serviceRegistryActor")
      val supervisor = TestProbe("supervisor")
      val copyActor = buildCopyActor(
        workflow = workflow,
        blacklistCache = blacklistCache,
        fakeIoActor = ioActor,
        fakeServiceRegistryActor = serviceRegistryActor,
        supervisor = supervisor,
        grouping = grouping
      )

      supervisor watch copyActor

      val command = buildCopyCommand(hitId = 3, bucket = WideOpenBucket)
      copyActor ! command

      eventually {
        copyActor.underlyingActor.stateName shouldBe WaitingForIoResponses
      }

      ioActor.expectMsgPF(5 seconds) { case ioCommand: DefaultIoCopyCommand =>
        ioActor.reply(IoFailure(ioCommand, new RuntimeException()))
      }

      val List(readHit, readBucket, writeHit) =
        instrumentationCounts(n = 3, serviceRegistryActor = serviceRegistryActor)

      readHit.bucket.path.toList shouldBe expectedMetric(Hit, Read, UntestedCacheResult)
      readBucket.bucket.path.toList shouldBe expectedMetric(Bucket, Read, GoodCacheResult)
      writeHit.bucket.path.toList shouldBe expectedMetric(Hit, Write, BadCacheResult)

      // Assert blacklist entries were made for bucket and hit.
      blacklistCache.bucketCache.size() shouldBe 2
      // The hit is bad but because the failure was generic the bucket which was marked good should stay that way.
      blacklistCache.bucketCache.get(WideOpenBucket) shouldBe GoodCacheResult
      blacklistCache.bucketCache.get(LockedDownBucket) shouldBe BadCacheResult

      blacklistCache.hitCache.size() shouldBe 4
      blacklistCache.hitCache.get(CallCachingEntryId(0)) shouldBe GoodCacheResult
      blacklistCache.hitCache.get(CallCachingEntryId(1)) shouldBe BadCacheResult
      blacklistCache.hitCache.get(CallCachingEntryId(2)) shouldBe UntestedCacheResult
      blacklistCache.hitCache.get(CallCachingEntryId(3)) shouldBe BadCacheResult

      supervisor.expectMsgPF(5 seconds) { case _: CopyingOutputsFailedResponse => }

      supervisor.expectTerminated(copyActor, 5 seconds)
    }

    val workflow2 = buildWorkflow(grouping)

    {
      // Step 4: a new workflow from the same grouping tries to copy the bad cache hit from step 3.
      val ioActor = TestProbe("ioActor")
      val serviceRegistryActor = TestProbe("serviceRegistryActor")
      val supervisor = TestProbe("supervisor")
      val copyActor = buildCopyActor(
        workflow = workflow2,
        blacklistCache = blacklistCache,
        fakeIoActor = ioActor,
        fakeServiceRegistryActor = serviceRegistryActor,
        supervisor = supervisor,
        grouping = grouping
      )

      supervisor watch copyActor

      val command = buildCopyCommand(hitId = 3, bucket = WideOpenBucket)
      copyActor ! command

      supervisor.expectMsgPF(5 seconds) { case _: CopyingOutputsFailedResponse =>
      }
      // The IoActor should not be consulted and the copying actor should simply stop itself without transitioning.
      supervisor.expectTerminated(copyActor)
      ioActor.expectNoMessage(5 seconds)

      val List(readHit) = instrumentationCounts(n = 1, serviceRegistryActor = serviceRegistryActor)

      readHit.bucket.path.toList shouldBe expectedMetric(Hit, Read, BadCacheResult)

      // Assert blacklist entries were made for bucket and hit.
      blacklistCache.bucketCache.size() shouldBe 2
      // The hit is bad but because the failure was generic the bucket which was previously marked good should stay that way.
      blacklistCache.bucketCache.get(WideOpenBucket) shouldBe GoodCacheResult
      blacklistCache.bucketCache.get(LockedDownBucket) shouldBe BadCacheResult

      blacklistCache.hitCache.size() shouldBe 4
      blacklistCache.hitCache.get(CallCachingEntryId(0)) shouldBe GoodCacheResult
      blacklistCache.hitCache.get(CallCachingEntryId(1)) shouldBe BadCacheResult
      blacklistCache.hitCache.get(CallCachingEntryId(2)) shouldBe UntestedCacheResult
      blacklistCache.hitCache.get(CallCachingEntryId(3)) shouldBe BadCacheResult
    }

    {
      // Step 5: a new workflow from the same grouping tries to copy an unknown cache hit from the bucket blacklisted in step 1.
      val ioActor = TestProbe("ioActor")
      val serviceRegistryActor = TestProbe("serviceRegistryActor")
      val supervisor = TestProbe("supervisor")
      val copyActor = buildCopyActor(
        workflow = workflow2,
        blacklistCache = blacklistCache,
        fakeIoActor = ioActor,
        fakeServiceRegistryActor = serviceRegistryActor,
        supervisor = supervisor,
        grouping = grouping
      )

      supervisor watch copyActor

      val command = buildCopyCommand(hitId = 4, bucket = LockedDownBucket)
      copyActor ! command

      supervisor.expectMsgPF(5 seconds) { case _: CopyingOutputsFailedResponse =>
      }
      // The IoActor should not be consulted and the copying actor should simply stop itself without transitioning.
      supervisor.expectTerminated(copyActor)
      ioActor.expectNoMessage(5 seconds)

      val List(readHit, readBucket) = instrumentationCounts(n = 2, serviceRegistryActor = serviceRegistryActor)

      readHit.bucket.path.toList shouldBe expectedMetric(Hit, Read, UntestedCacheResult)
      readBucket.bucket.path.toList shouldBe expectedMetric(Bucket, Read, BadCacheResult)

      // Assert blacklist entries were made for bucket and hit.
      blacklistCache.bucketCache.size() shouldBe 2
      // The hit is bad but because the failure was generic the bucket which was previously marked good should stay that way.
      blacklistCache.bucketCache.get(WideOpenBucket) shouldBe GoodCacheResult
      blacklistCache.bucketCache.get(LockedDownBucket) shouldBe BadCacheResult

      blacklistCache.hitCache.size() shouldBe 5
      blacklistCache.hitCache.get(CallCachingEntryId(0)) shouldBe GoodCacheResult
      blacklistCache.hitCache.get(CallCachingEntryId(1)) shouldBe BadCacheResult
      blacklistCache.hitCache.get(CallCachingEntryId(2)) shouldBe UntestedCacheResult
      blacklistCache.hitCache.get(CallCachingEntryId(3)) shouldBe BadCacheResult
      blacklistCache.hitCache.get(CallCachingEntryId(4)) shouldBe UntestedCacheResult
    }
  }

  private def instrumentationCounts(n: Int, serviceRegistryActor: TestProbe): List[CromwellCount] = {
    val received = serviceRegistryActor.receiveN(n = n, max = 5 seconds).toList
    val instrumentationCounts = received collect { case InstrumentationServiceMessage(c) => c } collect {
      case c: CromwellCount => c
    }
    instrumentationCounts foreach { c => c.value shouldBe 1; c.sampling shouldBe 1.0 }

    instrumentationCounts
  }

  type TestFSMRefBatchBackendCacheHitCopyingActor = TestFSMRef[
    StandardCacheHitCopyingActorState,
    Option[StandardCacheHitCopyingActorData],
    BatchBackendCacheHitCopyingActor
  ]

  private def buildWorkflow(grouping: Option[String]): HasWorkflowIdAndSources = {
    val workflowId = WorkflowId.randomId()
    new HasWorkflowIdAndSources {
      override val sources: WorkflowSourceFilesCollection = {
        val collection = mock[WorkflowSourceFilesCollection]
        val workflowOptions = grouping match {
          case None => WorkflowOptions.empty
          case Some(g) => WorkflowOptions.fromMap(Map("google_project" -> g)).get
        }
        collection.workflowOptions returns workflowOptions

        collection
      }

      override def id: WorkflowId = workflowId
    }
  }

  private def buildCopyActor(workflow: HasWorkflowIdAndSources,
                             blacklistCache: BlacklistCache,
                             fakeIoActor: TestProbe,
                             fakeServiceRegistryActor: TestProbe,
                             supervisor: TestProbe,
                             grouping: Option[String]
  ): TestFSMRefBatchBackendCacheHitCopyingActor = {
    // Couldn't mock this, possibly due to the use of `Refined` in two parameters:
    //
    // Underlying exception : java.lang.IllegalArgumentException: Cannot cast to primitive type: int
    // org.mockito.exceptions.base.MockitoException:
    // Mockito cannot mock this class: class cromwell.backend.google.batch.actors.GcpBatchConfigurationAttributes.
    val batchConfigurationAttributes = GcpBatchConfigurationAttributes(
      project = null,
      computeServiceAccount = null,
      auths = null,
      dockerhubToken = null,
      enableFuse = false,
      executionBucket = null,
      location = null,
      maxPollingInterval = 0,
      qps = refineMV[Positive](10),
      cacheHitDuplicationStrategy = CopyCachedOutputs,
      requestWorkers = refineMV[Positive](1),
      batchTimeout = null,
      dockerMirroringOpt = None,
      logFlushPeriod = None,
      gcsTransferConfiguration = null,
      virtualPrivateCloudConfiguration = VirtualPrivateCloudConfiguration(None, None),
      batchRequestTimeoutConfiguration = null,
      referenceFileToDiskImageMappingOpt = None,
      checkpointingInterval = 10.minutes,
      logsPolicy = GcpBatchLogsPolicy.CloudLogging,
      maxTransientErrorRetries = 10
    )

    val batchConfiguration = mockWithDefaults[GcpBatchConfiguration]
    batchConfiguration.batchAttributes returns batchConfigurationAttributes

    val commandTaskDefinition = mock[CommandTaskDefinition]
    commandTaskDefinition.outputs returns List.empty
    val commandCallNode = CommandCallNode(
      identifier = WomIdentifier(LocalName("bar"), FullyQualifiedName("foo.bar")),
      callable = commandTaskDefinition,
      inputPorts = Set.empty,
      inputDefinitionMappings = List.empty,
      nonInputBasedPrerequisites = Set.empty,
      outputIdentifierCompoundingFunction = null,
      sourceLocation = None
    )

    val backendJobDescriptorKey = BackendJobDescriptorKey(
      call = commandCallNode,
      index = None,
      attempt = 1
    )

    // noinspection ScalaUnusedSymbol
    def mapper(jobPaths: GcpBatchJobPaths, originalPath: String): String = originalPath

    val workflowDescriptor = mock[BackendWorkflowDescriptor]
    workflowDescriptor.id returns workflow.id
    workflowDescriptor.workflowOptions returns WorkflowOptions.fromMap(Map("batch_gcs_root" -> "foo")).get

    val workflowPaths = mock[GcpBatchWorkflowPaths]
    workflowPaths.standardStreamNameToFileNameMetadataMapper returns mapper
    workflowPaths.workflowRoot returns mock[Path]
    val batchJobPaths = mock[GcpBatchJobPaths]
    batchJobPaths.workflowPaths returns workflowPaths

    val copyDestinationPaths = mock[GcpBatchJobPaths]
    val copyDestinationRcPath = mock[Path]
    copyDestinationPaths.detritusPaths returns Map(JobPaths.ReturnCodePathKey -> copyDestinationRcPath)

    batchJobPaths.forCallCacheCopyAttempts returns copyDestinationPaths
    batchJobPaths.metadataPaths returns Map.empty
    workflowPaths.toJobPaths(any[BackendJobDescriptor]).returns(batchJobPaths)

    def identityPathMocker(str: Any): Try[Path] = {
      val path = mock[Path]
      path.toString returns str.asInstanceOf[String]
      Success(path)
    }

    workflowPaths.getPath(anyString).answers(identityPathMocker)

    val runtimeAttributesBuilder = mock[StandardValidatedRuntimeAttributesBuilder]
    runtimeAttributesBuilder
      .build(any[Map[String, WomValue]], any[Logger])
      .returns(ValidatedRuntimeAttributes(Map.empty))

    val backendInitializationData = mock[GcpBackendInitializationData]
    backendInitializationData.gcpBatchConfiguration returns batchConfiguration
    backendInitializationData.workflowPaths returns workflowPaths
    backendInitializationData.runtimeAttributesBuilder returns runtimeAttributesBuilder

    val backendJobDescriptor = BackendJobDescriptor(
      workflowDescriptor = workflowDescriptor,
      key = backendJobDescriptorKey,
      runtimeAttributes = Map.empty,
      evaluatedTaskInputs = Map.empty,
      maybeCallCachingEligible = DockerWithHash("foo"),
      dockerSize = None,
      prefetchedKvStoreEntries = Map.empty
    )

    val params = DefaultStandardCacheHitCopyingActorParams(
      jobDescriptor = backendJobDescriptor,
      backendInitializationDataOption = Option(backendInitializationData),
      serviceRegistryActor = fakeServiceRegistryActor.ref,
      ioActor = fakeIoActor.ref,
      configurationDescriptor = null,
      cacheCopyAttempt = 0,
      blacklistCache = Option(blacklistCache)
    )

    val actorUnderTest = TestFSMRef(new BatchBackendCacheHitCopyingActor(params), supervisor = supervisor.ref)

    eventually {
      actorUnderTest.underlyingActor.stateName shouldBe Idle
    }

    actorUnderTest
  }

  private def buildCopyCommand(hitId: Long, bucket: String): CopyOutputsCommand = {
    val callRoot = s"gs://$bucket/workflow-id/call-name"
    val rcFile = callRoot + "/rc"

    CopyOutputsCommand(
      womValueSimpletons = List.empty,
      jobDetritusFiles = Map(JobPaths.CallRootPathKey -> callRoot, JobPaths.ReturnCodePathKey -> rcFile),
      returnCode = Option(0),
      cacheHit = CallCachingEntryId(hitId)
    )
  }

  sealed trait BlacklistingType extends HasFormatting
  case object Hit extends BlacklistingType
  case object Bucket extends BlacklistingType

  sealed trait CacheAccessType extends HasFormatting
  case object Read extends CacheAccessType
  case object Write extends CacheAccessType

  private def expectedMetric(hitOrBucket: BlacklistingType,
                             accessType: CacheAccessType,
                             status: BlacklistStatus
  ): List[String] =
    List("job",
         "callcaching",
         "blacklist",
         accessType.metricFormat,
         hitOrBucket.metricFormat,
         status.getClass.getSimpleName.dropRight(1)
    )
}
