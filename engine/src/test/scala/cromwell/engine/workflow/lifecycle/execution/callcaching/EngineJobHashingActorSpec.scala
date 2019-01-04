package cromwell.engine.workflow.lifecycle.execution.callcaching

import wdl.draft2.model.command.StringCommandPart
import akka.actor.{Actor, ActorRef, Props}
import akka.testkit.{TestActorRef, TestProbe}
import cats.syntax.validated._
import cromwell.backend._
import cromwell.core._
import cromwell.core.callcaching._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheHashingJobActor.{CompleteFileHashingResult, InitialHashingResult, NoFileHashesResult}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadingJobActor.NextHit
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor._
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.util.WomMocks
import org.scalatest.concurrent.Eventually
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpecLike, Matchers}
import wom.core.LocallyQualifiedName
import wom.graph.WomIdentifier
import wom.values.WomValue

class EngineJobHashingActorSpec extends TestKitSuite with FlatSpecLike with Matchers with BackendSpec with TableDrivenPropertyChecks with Eventually {
  behavior of "EngineJobHashingActor"

  def templateJobDescriptor(inputs: Map[LocallyQualifiedName, WomValue] = Map.empty) = {
    val task = WomMocks.mockTaskDefinition("hello").copy(
      commandTemplateBuilder = Function.const(List(StringCommandPart("Do the stuff... now!!")).validNel)
    )
    val call = WomMocks.mockTaskCall(WomIdentifier("hello", "workflow.hello")).copy(callable = task)
    val workflowDescriptor = mock[BackendWorkflowDescriptor]
    workflowDescriptor.id returns WorkflowId.randomId()
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, BackendJobDescriptorKey(call, None, 1), Map.empty, fqnWdlMapToDeclarationMap(inputs), NoDocker, None, Map.empty)
    jobDescriptor
  }
  
  val serviceRegistryActorProbe = TestProbe()

  def makeEJHA(receiver: ActorRef, activity: CallCachingActivity, ccReaderProps: Props = Props.empty) = {
    TestActorRef[EngineJobHashingActor](
      EngineJobHashingActorTest.props(
        receiver,
        serviceRegistryActorProbe.ref,
        templateJobDescriptor(),
        None,
        Props.empty,
        ccReaderProps,
        Set.empty,
        "backend",
        activity,
        DockerWithHash("ubuntu@sha256:blablabla")
      )
    )
  }

  it should "record initial hashes" in {
    val receiver = TestProbe()
    val activity = CallCachingActivity(ReadAndWriteCache)
    val actorUnderTest = makeEJHA(receiver.ref, activity)

    val initialResult: InitialHashingResult = mock[InitialHashingResult]
    actorUnderTest ! initialResult
    eventually {
      actorUnderTest.underlyingActor.initialHash shouldBe Some(initialResult)
    }
  }

  it should "create a CCReader actor or not depending on CC activity" in {
    val activities = Table(
      ("readWriteMode", "hasCCReadActor"),
      (ReadCache, true),
      (WriteCache, false),
      (ReadAndWriteCache, true)
    )
    forAll(activities) { case (readWriteMode, hasCCReadActor) =>
      val receiver = TestProbe()
      val actorUnderTest = makeEJHA(receiver.ref, CallCachingActivity(readWriteMode))
      actorUnderTest.underlyingActor.callCacheReadingJobActor.isDefined shouldBe hasCCReadActor
    }
  }

  it should "send hashes to receiver when receiving a NoFileHashesResult" in {
    val receiver = TestProbe()
    val actorUnderTest = makeEJHA(receiver.ref, CallCachingActivity(ReadAndWriteCache))
    val initialHashes = Set(HashResult(HashKey("key"), HashValue("value")))
    val initialAggregatedHash = "aggregatedHash"
    val initialResult = InitialHashingResult(initialHashes, initialAggregatedHash)
    actorUnderTest ! initialResult
    actorUnderTest ! NoFileHashesResult
    receiver.expectMsg(CallCacheHashes(initialHashes, initialAggregatedHash, None))
  }

  it should "send hashes to receiver when receiving a CompleteFileHashingResult" in {
    val receiver = TestProbe()
    val actorUnderTest = makeEJHA(receiver.ref, CallCachingActivity(ReadAndWriteCache))
    
    val initialHashes = Set(HashResult(HashKey("key"), HashValue("value")))
    val initialAggregatedHash = "aggregatedHash"
    val initialResult = InitialHashingResult(initialHashes, initialAggregatedHash)
    val fileHashes = Set(HashResult(HashKey("file key"), HashValue("value")))
    val fileAggregatedHash = "aggregatedFileHash"
    val fileResult = CompleteFileHashingResult(fileHashes, fileAggregatedHash)
    
    actorUnderTest ! initialResult
    actorUnderTest ! fileResult
    receiver.expectMsg(CallCacheHashes(initialHashes, initialAggregatedHash, Option(FileHashes(fileHashes, fileAggregatedHash))))
  }

  it should "forward CacheMiss to receiver" in {
    val receiver = TestProbe()
    val activity = CallCachingActivity(ReadAndWriteCache)
    val actorUnderTest = makeEJHA(receiver.ref, activity)

    actorUnderTest ! CacheMiss
    receiver.expectMsg(CacheMiss)
  }

  it should "forward CacheHit to receiver" in {
    val receiver = TestProbe()
    val activity = CallCachingActivity(ReadAndWriteCache)
    val actorUnderTest = makeEJHA(receiver.ref, activity)

    val cacheHit = mock[CacheHit]
    actorUnderTest ! cacheHit
    receiver.expectMsg(cacheHit)
  }

  it should "forward NextHit to CCRead actor" in {
    val receiver = TestProbe()
    val activity = CallCachingActivity(ReadAndWriteCache)
    val monitorProbe = TestProbe()
    val ccReadActorProps = Props(new Actor {
      override def receive: Receive = {
        case NextHit => monitorProbe.ref forward NextHit
      }
    })
    
    val actorUnderTest = makeEJHA(receiver.ref, activity, ccReadActorProps)

    actorUnderTest ! NextHit
    monitorProbe.expectMsg(NextHit)
  }

  it should "fail if it receives NextHit and doesn't have a CCRead actor" in {
    val receiver = TestProbe()
    val activity = CallCachingActivity(WriteCache)
    val actorUnderTest = makeEJHA(receiver.ref, activity)
    receiver.watch(actorUnderTest)
    actorUnderTest ! NextHit
    receiver.expectMsgClass(classOf[HashError])
    serviceRegistryActorProbe.expectMsgClass(classOf[PutMetadataAction])
    receiver.expectTerminated(actorUnderTest)
  }

  it should "fail if it receives a HashingFailedMessage" in {
    val receiver = TestProbe()
    val activity = CallCachingActivity(ReadAndWriteCache)
    val actorUnderTest = makeEJHA(receiver.ref, activity)
    receiver.watch(actorUnderTest)
    actorUnderTest ! HashingFailedMessage("someFile", new Exception("[TEST] Some exception"))
    receiver.expectMsgClass(classOf[HashError])
    serviceRegistryActorProbe.expectMsgClass(classOf[PutMetadataAction])
    receiver.expectTerminated(actorUnderTest)
  }

  it should "fail if it receives a FinalFileHashingResult but has no InitialHashingResult" in {
    val receiver = TestProbe()
    val activity = CallCachingActivity(ReadAndWriteCache)
    val actorUnderTest = makeEJHA(receiver.ref, activity)
    receiver.watch(actorUnderTest)
    actorUnderTest ! NoFileHashesResult
    receiver.expectMsgClass(classOf[HashError])
    serviceRegistryActorProbe.expectMsgClass(classOf[PutMetadataAction])
    receiver.expectTerminated(actorUnderTest)
  }
  
  object EngineJobHashingActorTest {
    def props(receiver: ActorRef,
              serviceRegistryActor: ActorRef,
              jobDescriptor: BackendJobDescriptor,
              initializationData: Option[BackendInitializationData],
              fileHashingActorProps: Props,
              callCacheReadingJobActorProps: Props,
              runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition],
              backendName: String,
              activity: CallCachingActivity,
              callCachingEligible: CallCachingEligible): Props = Props(new EngineJobHashingActorTest(
      receiver = receiver,
      serviceRegistryActor = serviceRegistryActor,
      jobDescriptor = jobDescriptor,
      initializationData = initializationData,
      fileHashingActorProps = fileHashingActorProps,
      callCacheReadingJobActorProps = callCacheReadingJobActorProps,
      runtimeAttributeDefinitions = runtimeAttributeDefinitions,
      backendName = backendName,
      activity = activity,
      callCachingEligible = callCachingEligible))
  }
  
  class EngineJobHashingActorTest(receiver: ActorRef,
                                  serviceRegistryActor: ActorRef,
                                  jobDescriptor: BackendJobDescriptor,
                                  initializationData: Option[BackendInitializationData],
                                  fileHashingActorProps: Props,
                                  callCacheReadingJobActorProps: Props,
                                  runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition],
                                  backendName: String,
                                  activity: CallCachingActivity,
                                  callCachingEligible: CallCachingEligible) extends EngineJobHashingActor(
    receiver = receiver,
    serviceRegistryActor = serviceRegistryActor,
    jobDescriptor = jobDescriptor,
    initializationData = initializationData,
    fileHashingActorProps = fileHashingActorProps,
    callCacheReadingJobActorProps = callCacheReadingJobActorProps,
    runtimeAttributeDefinitions = runtimeAttributeDefinitions,
    backendNameForCallCachingPurposes = backendName,
    activity = activity,
    callCachingEligible = callCachingEligible,
    callCachePathPrefixes = None) {
    // override preStart to nothing to prevent the creation of the CCHJA.
    // This way it doesn't interfere with the tests and we can manually inject the messages we want
    override def preStart() =  ()
  }

}
