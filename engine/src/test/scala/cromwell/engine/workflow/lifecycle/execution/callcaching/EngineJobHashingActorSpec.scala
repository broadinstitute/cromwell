package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{ActorRef, ActorSystem, Props}
import akka.testkit.{ImplicitSender, TestKit, TestProbe}
import cromwell.CromwellTestkitSpec
import cromwell.backend.callcaching.FileHasherWorkerActor.{FileHashResponse, SingleFileHashRequest}
import cromwell.backend.{BackendInitializationData, BackendJobDescriptor, BackendJobDescriptorKey, BackendWorkflowDescriptor, RuntimeAttributeDefinition}
import cromwell.core.callcaching._
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.{CacheHit, CacheMiss, CallCacheHashes}
import org.scalatest.mockito.MockitoSugar
import org.scalatest.{BeforeAndAfterAll, Matchers, WordSpecLike}
import wdl4s._
import wdl4s.command.CommandPart
import wdl4s.values.{WdlFile, WdlValue}

import scala.concurrent.duration._
import scala.language.postfixOps

class EngineJobHashingActorSpec extends TestKit(new CromwellTestkitSpec.TestWorkflowManagerSystem().actorSystem)
  with ImplicitSender with WordSpecLike with Matchers with MockitoSugar with BeforeAndAfterAll {

  import EngineJobHashingActorSpec._

  implicit val actorSystem: ActorSystem = system

  val readModes = List(CallCachingActivity(ReadCache), CallCachingActivity(ReadAndWriteCache))
  val writeModes = List(CallCachingActivity(WriteCache), CallCachingActivity(ReadAndWriteCache))
  val allModes = List(CallCachingActivity(ReadCache), CallCachingActivity(WriteCache), CallCachingActivity(ReadAndWriteCache))

  "Engine job hashing actor" must {
    allModes foreach { activity =>
      val expectation = activity.readWriteMode match {
        case ReadCache => "cache hit"
        case WriteCache => "hashes"
        case ReadAndWriteCache => "cache hit and hashes"
      }

      s"Respect the CallCachingMode and report back $expectation for the ${activity.readWriteMode} activity" in {
        val singleMetaInfoIdSet = Set(MetaInfoId(1))
        val replyTo = TestProbe()
        val deathWatch = TestProbe()

        val cacheLookupResponses: Map[String, Set[MetaInfoId]] = if (activity.readFromCache) standardCacheLookupResponses(singleMetaInfoIdSet, singleMetaInfoIdSet, singleMetaInfoIdSet) else Map.empty
        val ejha = createEngineJobHashingActor(
          replyTo = replyTo.ref,
          activity = activity,
          cacheLookupResponses = cacheLookupResponses)

        deathWatch watch ejha

        if (activity.readFromCache) replyTo.expectMsg(CacheHit(MetaInfoId(1)))
        if (activity.writeToCache) replyTo.expectMsgPF(max = 5 seconds, hint = "awaiting cache hit message") {
          case CallCacheHashes(hashes) => hashes.size should be(3)
          case x => fail(s"Cache hit anticipated! Instead got a ${x.getClass.getSimpleName}")
        }

        deathWatch.expectTerminated(ejha, 5 seconds)
      }

      s"Wait for requests to the FileHasherActor for the ${activity.readWriteMode} activity" in {
        val singleMetaInfoIdSet = Set(MetaInfoId(1))
        val replyTo = TestProbe()
        val fileHasherActor = TestProbe()
        val deathWatch = TestProbe()

        val initialCacheLookupResponses: Map[String, Set[MetaInfoId]] = if (activity.readFromCache) standardCacheLookupResponses(singleMetaInfoIdSet, singleMetaInfoIdSet, singleMetaInfoIdSet) else Map.empty
        val fileCacheLookupResponses = Map("input: File inputFile1" -> singleMetaInfoIdSet, "input: File inputFile2" -> singleMetaInfoIdSet)

        val jobDescriptor = templateJobDescriptor(inputs = Map(
            "inputFile1" -> WdlFile("path"),
            "inputFile2" -> WdlFile("path")))
        val ejha = createEngineJobHashingActor(
          replyTo = replyTo.ref,
          activity = activity,
          jobDescriptor = jobDescriptor,
          fileHashingActor = Some(fileHasherActor.ref),
          cacheLookupResponses = initialCacheLookupResponses ++ fileCacheLookupResponses)

        deathWatch watch ejha

        twice { iteration =>
          fileHasherActor.expectMsgPF(max = 5 seconds, hint = s"awaiting file hash request #$iteration") {
            case SingleFileHashRequest(jobKey, hashKey, file, initializationData) =>
              file should be(WdlFile("path"))
              fileHasherActor.send(ejha, FileHashResponse(HashResult(hashKey, HashValue("blah di blah"))))
            case x => fail(s"SingleFileHashRequest anticipated! Instead got a ${x.getClass.getSimpleName}")
          }
        }

        if (activity.readFromCache) replyTo.expectMsg(CacheHit(MetaInfoId(1)))
        if (activity.writeToCache) replyTo.expectMsgPF(max = 5 seconds, hint = "awaiting cache hit message") {
          case CallCacheHashes(hashes) => hashes.size should be(5)
          case x => fail(s"Cache hit anticipated! Instead got a ${x.getClass.getSimpleName}")
        }

        deathWatch.expectTerminated(ejha, 5 seconds)
      }

      s"Cache miss for bad FileHasherActor results but still return hashes in the ${activity.readWriteMode} activity" in {
        val singleMetaInfoIdSet = Set(MetaInfoId(1))
        val replyTo = TestProbe()
        val fileHasherActor = TestProbe()
        val deathWatch = TestProbe()

        val initialCacheLookupResponses: Map[String, Set[MetaInfoId]] = if (activity.readFromCache) standardCacheLookupResponses(singleMetaInfoIdSet, singleMetaInfoIdSet, singleMetaInfoIdSet) else Map.empty
        val fileCacheLookupResponses = Map("input: File inputFile1" -> Set(MetaInfoId(2)), "input: File inputFile2" -> singleMetaInfoIdSet)

        val jobDescriptor = templateJobDescriptor(inputs = Map(
          "inputFile1" -> WdlFile("path"),
          "inputFile2" -> WdlFile("path")))
        val ejha = createEngineJobHashingActor(
          replyTo = replyTo.ref,
          activity = activity,
          jobDescriptor = jobDescriptor,
          fileHashingActor = Some(fileHasherActor.ref),
          cacheLookupResponses = initialCacheLookupResponses ++ fileCacheLookupResponses)

        deathWatch watch ejha

        // Hello, future Cromwellian! I imagine you're reading this because you've just introduced file hash short-circuiting on cache miss and,
        // depending on timings, you may not get the second file hash request in read-only mode. You might want to refactor this test to reply
        // only to the "cache miss" file and then check that the test probe receives the appropriate "cancellation" message.
        // ... or not! Don't just blindly allow a ghost of the past to tell you what to do! Live your own life and excel!
        twice { iteration =>
          fileHasherActor.expectMsgPF(max = 5 seconds, hint = s"awaiting file hash request #$iteration") {
            case SingleFileHashRequest(jobKey, hashKey, file, initializationData) =>
              file should be(WdlFile("path"))
              fileHasherActor.send(ejha, FileHashResponse(HashResult(hashKey, HashValue("blah di blah"))))
            case x => fail(s"SingleFileHashRequest anticipated! Instead got a ${x.getClass.getSimpleName}")
          }
        }

        if (activity.readFromCache) replyTo.expectMsg(CacheMiss)
        if (activity.writeToCache) replyTo.expectMsgPF(max = 5 seconds, hint = "awaiting cache hit message") {
          case CallCacheHashes(hashes) => hashes.size should be(5)
          case x => fail(s"Cache hit anticipated! Instead got a ${x.getClass.getSimpleName}")
        }

        deathWatch.expectTerminated(ejha, 5 seconds)
      }

      s"Detect call cache misses for the ${activity.readWriteMode} activity" in {
        val singleMetaInfoIdSet = Set(MetaInfoId(1))
        val replyTo = TestProbe()
        val deathWatch = TestProbe()

        val cacheLookupResponses: Map[String, Set[MetaInfoId]] = if (activity.readFromCache) standardCacheLookupResponses(singleMetaInfoIdSet, singleMetaInfoIdSet, Set(MetaInfoId(2))) else Map.empty
        val ejha = createEngineJobHashingActor(
          replyTo = replyTo.ref,
          activity = activity,
          cacheLookupResponses = cacheLookupResponses)

        deathWatch watch ejha

        if (activity.readFromCache) replyTo.expectMsg(CacheMiss)
        if (activity.writeToCache) replyTo.expectMsgPF(max = 5 seconds, hint = "awaiting cache hit message") {
          case CallCacheHashes(hashes) => hashes.size should be(3)
          case x => fail(s"Cache hit anticipated! Instead got a ${x.getClass.getSimpleName}")
        }

        deathWatch.expectTerminated(ejha, 5 seconds)
      }
    }
  }

  override def afterAll() = {
    TestKit.shutdownActorSystem(system)
  }
}

object EngineJobHashingActorSpec extends MockitoSugar {
  import org.mockito.Mockito._

  def createEngineJobHashingActor
  (
    replyTo: ActorRef,
    activity: CallCachingActivity,
    jobDescriptor: BackendJobDescriptor = templateJobDescriptor(),
    initializationData: Option[BackendInitializationData] = None,
    fileHashingActor: Option[ActorRef] = None,
    cacheLookupResponses: Map[String, Set[MetaInfoId]] = Map.empty,
    runtimeAttributeDefinitions: Set[RuntimeAttributeDefinition] = Set.empty,
    backendName: String = "whatever"
  )(implicit system: ActorSystem) = {
    val callCacheReadActor = system.actorOf(Props(new PredictableCallCacheReadActor(cacheLookupResponses)))
    system.actorOf(EngineJobHashingActor.props(
      receiver = replyTo,
      jobDescriptor = jobDescriptor,
      initializationData = initializationData,
      fileHashingActor = fileHashingActor.getOrElse(emptyActor),
      callCacheReadActor = callCacheReadActor,
      runtimeAttributeDefinitions = runtimeAttributeDefinitions,
      backendName = backendName,
      activity = activity))
  }

  def emptyActor(implicit actorSystem: ActorSystem) = actorSystem.actorOf(Props.empty)

  def templateJobDescriptor(inputs: Map[LocallyQualifiedName, WdlValue] = Map.empty) = {
    val commandPart = mock[CommandPart]
    val task = mock[Task]
    val call = mock[Call]
    when(task.commandTemplate).thenReturn(Seq(commandPart))
    when(task.outputs).thenReturn(List.empty)
    when(call.task).thenReturn(task)
    when(commandPart.toString).thenReturn("literally don't care")
    val workflowDescriptor = mock[BackendWorkflowDescriptor]
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, BackendJobDescriptorKey(call, None, 1), Map.empty, inputs)
    jobDescriptor
  }

  def standardCacheLookupResponses(commandTemplate: Set[MetaInfoId], inputCount: Set[MetaInfoId], backendName: Set[MetaInfoId]) = Map(
    "command template" -> commandTemplate,
    "input count" -> inputCount,
    "backend name" -> backendName
  )

  def twice[A](block: Int => A) = {
    block(1)
    block(2)
  }
}