package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{Actor, ActorRef, Props}
import akka.testkit.{TestActorRef, TestProbe}
import cromwell.backend._
import cromwell.core._
import cromwell.core.callcaching._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheHashingJobActor.{CompleteFileHashingResult, InitialHashingResult, NoFileHashesResult}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadingJobActor.NextHit
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor._
import org.scalatest.concurrent.Eventually
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpecLike, Matchers}
import wdl4s.values.WdlValue
import wdl4s.{Task, TaskCall}

class EngineJobHashingActorSpec extends TestKitSuite with FlatSpecLike with Matchers with BackendSpec with TableDrivenPropertyChecks with Eventually {
  behavior of "EngineJobHashingActor"

  def templateJobDescriptor(inputs: Map[LocallyQualifiedName, WdlValue] = Map.empty) = {
    val task = mock[Task]
    val call = mock[TaskCall]
    task.commandTemplateString returns "Do the stuff... now!!"
    task.outputs returns List.empty
    task.fullyQualifiedName returns "workflow.hello"
    call.task returns task
    val workflowDescriptor = mock[BackendWorkflowDescriptor]
    workflowDescriptor.id returns WorkflowId.randomId()
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, BackendJobDescriptorKey(call, None, 1), Map.empty, fqnMapToDeclarationMap(inputs), NoDocker, Map.empty)
    jobDescriptor
  }

  def makeEJHA(receiver: ActorRef, activity: CallCachingActivity, ccReaderProps: Props = Props.empty) = {
    TestActorRef[EngineJobHashingActor](
      EngineJobHashingActorTest.props(
        receiver,
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
    forAll(activities) { case ((readWriteMode, hasCCReadActor)) =>
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
    receiver.expectTerminated(actorUnderTest)
  }

  it should "fail if it receives a HashingFailedMessage" in {
    val receiver = TestProbe()
    val activity = CallCachingActivity(ReadAndWriteCache)
    val actorUnderTest = makeEJHA(receiver.ref, activity)
    receiver.watch(actorUnderTest)
    actorUnderTest ! mock[HashingFailedMessage]
    receiver.expectMsgClass(classOf[HashError])
    receiver.expectTerminated(actorUnderTest)
  }

  it should "fail if it receives a FinalFileHashingResult but has no InitialHashingResult" in {
    val receiver = TestProbe()
    val activity = CallCachingActivity(ReadAndWriteCache)
    val actorUnderTest = makeEJHA(receiver.ref, activity)
    receiver.watch(actorUnderTest)
    actorUnderTest ! NoFileHashesResult
    receiver.expectMsgClass(classOf[HashError])
    receiver.expectTerminated(actorUnderTest)
  }
  
  object EngineJobHashingActorTest {
    def props(receiver: ActorRef,
              jobDescriptor: BackendJobDescriptor,
              initializationData: Option[BackendInitializationData],
              fileHashingActorProps: Props,
              callCacheReadingJobActorProps: Props,
              runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition],
              backendName: String,
              activity: CallCachingActivity,
              callCachingEligible: CallCachingEligible): Props = Props(new EngineJobHashingActorTest(
      receiver = receiver,
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
                                  jobDescriptor: BackendJobDescriptor,
                                  initializationData: Option[BackendInitializationData],
                                  fileHashingActorProps: Props,
                                  callCacheReadingJobActorProps: Props,
                                  runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition],
                                  backendName: String,
                                  activity: CallCachingActivity,
                                  callCachingEligible: CallCachingEligible) extends EngineJobHashingActor(
    receiver = receiver,
    jobDescriptor = jobDescriptor,
    initializationData = initializationData,
    fileHashingActorProps = fileHashingActorProps,
    callCacheReadingJobActorProps = callCacheReadingJobActorProps,
    runtimeAttributeDefinitions = runtimeAttributeDefinitions,
    backendName = backendName,
    activity = activity,
    callCachingEligible = callCachingEligible) {
    // override preStart to nothing to prevent the creation of the CCHJA.
    // This way it doesn't interfere with the tests and we can manually inject the messages we want
    override def preStart() =  ()
  }

}
