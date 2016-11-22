package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{ActorRef, ActorSystem, Props}
import akka.testkit.{ImplicitSender, TestProbe}
import cats.data.NonEmptyList
import cromwell.CromwellTestKitSpec
import cromwell.backend._
import cromwell.backend.callcaching.FileHashingActor.{FileHashResponse, SingleFileHashRequest}
import cromwell.core.callcaching._
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.{CacheHit, CacheMiss, CallCacheHashes}
import org.scalatest.mockito.MockitoSugar
import org.scalatest.{Matchers, WordSpecLike}
import wdl4s._
import wdl4s.values.{WdlFile, WdlValue}

import scala.concurrent.duration._
import scala.language.postfixOps

class EngineJobHashingActorSpec extends CromwellTestKitSpec
  with ImplicitSender with WordSpecLike with Matchers with MockitoSugar {

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
        val singleCallCachingEntryIdSet = Set(CallCachingEntryId(1))
        val replyTo = TestProbe()
        val deathWatch = TestProbe()

        val cacheLookupResponses: Map[String, Set[CallCachingEntryId]] = if (activity.readFromCache) standardCacheLookupResponses(singleCallCachingEntryIdSet, singleCallCachingEntryIdSet, singleCallCachingEntryIdSet, singleCallCachingEntryIdSet) else Map.empty
        val ejha = createEngineJobHashingActor(
          replyTo = replyTo.ref,
          activity = activity,
          cacheLookupResponses = cacheLookupResponses)

        deathWatch watch ejha

        if (activity.readFromCache) replyTo.expectMsg(CacheHit(NonEmptyList.of(CallCachingEntryId(1))))
        if (activity.writeToCache) replyTo.expectMsgPF(max = 5 seconds, hint = "awaiting cache hit message") {
          case CallCacheHashes(hashes) => hashes.size should be(4)
          case x => fail(s"Cache hit anticipated! Instead got a ${x.getClass.getSimpleName}")
        }

        deathWatch.expectTerminated(ejha, 5 seconds)
      }

      s"Wait for requests to the FileHashingActor for the ${activity.readWriteMode} activity" in {
        val singleCallCachingEntryIdSet = Set(CallCachingEntryId(1))
        val replyTo = TestProbe()
        val fileHashingActor = TestProbe()
        val deathWatch = TestProbe()

        val initialCacheLookupResponses: Map[String, Set[CallCachingEntryId]] = if (activity.readFromCache) standardCacheLookupResponses(singleCallCachingEntryIdSet, singleCallCachingEntryIdSet, singleCallCachingEntryIdSet, singleCallCachingEntryIdSet) else Map.empty
        val fileCacheLookupResponses = Map("input: File inputFile1" -> singleCallCachingEntryIdSet, "input: File inputFile2" -> singleCallCachingEntryIdSet)

        val jobDescriptor = templateJobDescriptor(inputs = Map(
            "inputFile1" -> WdlFile("path"),
            "inputFile2" -> WdlFile("path")))
        val ejha = createEngineJobHashingActor(
          replyTo = replyTo.ref,
          activity = activity,
          jobDescriptor = jobDescriptor,
          fileHashingActor = Option(fileHashingActor.ref),
          cacheLookupResponses = initialCacheLookupResponses ++ fileCacheLookupResponses)

        deathWatch watch ejha

        twice { iteration =>
          fileHashingActor.expectMsgPF(max = 5 seconds, hint = s"awaiting file hash request #$iteration") {
            case SingleFileHashRequest(jobKey, hashKey, file, initializationData) =>
              file should be(WdlFile("path"))
              fileHashingActor.send(ejha, FileHashResponse(HashResult(hashKey, HashValue("blah di blah"))))
            case x => fail(s"SingleFileHashRequest anticipated! Instead got a ${x.getClass.getSimpleName}")
          }
        }

        if (activity.readFromCache) replyTo.expectMsg(CacheHit(NonEmptyList.of(CallCachingEntryId(1))))
        if (activity.writeToCache) replyTo.expectMsgPF(max = 5 seconds, hint = "awaiting cache hit message") {
          case CallCacheHashes(hashes) => hashes.size should be(6)
          case x => fail(s"Cache hit anticipated! Instead got a ${x.getClass.getSimpleName}")
        }

        deathWatch.expectTerminated(ejha, 5 seconds)
      }

      s"Cache miss for bad FileHashingActor results but still return hashes in the ${activity.readWriteMode} activity" in {
        val singleCallCachingEntryIdSet = Set(CallCachingEntryId(1))
        val replyTo = TestProbe()
        val fileHashingActor = TestProbe()
        val deathWatch = TestProbe()

        val initialCacheLookupResponses: Map[String, Set[CallCachingEntryId]] = if (activity.readFromCache) standardCacheLookupResponses(singleCallCachingEntryIdSet, singleCallCachingEntryIdSet, singleCallCachingEntryIdSet, singleCallCachingEntryIdSet) else Map.empty
        val fileCacheLookupResponses = Map("input: File inputFile1" -> Set(CallCachingEntryId(2)), "input: File inputFile2" -> singleCallCachingEntryIdSet)

        val jobDescriptor = templateJobDescriptor(inputs = Map(
          "inputFile1" -> WdlFile("path"),
          "inputFile2" -> WdlFile("path")))
        val ejha = createEngineJobHashingActor(
          replyTo = replyTo.ref,
          activity = activity,
          jobDescriptor = jobDescriptor,
          fileHashingActor = Option(fileHashingActor.ref),
          cacheLookupResponses = initialCacheLookupResponses ++ fileCacheLookupResponses)

        deathWatch watch ejha

        // Hello, future Cromwellian! I imagine you're reading this because you've just introduced file hash short-circuiting on cache miss and,
        // depending on timings, you may not get the second file hash request in read-only mode. You might want to refactor this test to reply
        // only to the "cache miss" file and then check that the test probe receives the appropriate "cancellation" message.
        // ... or not! Don't just blindly allow a ghost of the past to tell you what to do! Live your own life and excel!
        twice { iteration =>
          fileHashingActor.expectMsgPF(max = 5 seconds, hint = s"awaiting file hash request #$iteration") {
            case SingleFileHashRequest(jobKey, hashKey, file, initializationData) =>
              file should be(WdlFile("path"))
              fileHashingActor.send(ejha, FileHashResponse(HashResult(hashKey, HashValue("blah di blah"))))
            case x => fail(s"SingleFileHashRequest anticipated! Instead got a ${x.getClass.getSimpleName}")
          }
        }

        if (activity.readFromCache) replyTo.expectMsg(CacheMiss)
        if (activity.writeToCache) replyTo.expectMsgPF(max = 5 seconds, hint = "awaiting cache hit message") {
          case CallCacheHashes(hashes) => hashes.size should be(6)
          case x => fail(s"Cache hit anticipated! Instead got a ${x.getClass.getSimpleName}")
        }

        deathWatch.expectTerminated(ejha, 5 seconds)
      }

      s"Detect call cache misses for the ${activity.readWriteMode} activity" in {
        val singleCallCachingEntryIdSet = Set(CallCachingEntryId(1))
        val replyTo = TestProbe()
        val deathWatch = TestProbe()

        val cacheLookupResponses: Map[String, Set[CallCachingEntryId]] = if (activity.readFromCache) standardCacheLookupResponses(singleCallCachingEntryIdSet, singleCallCachingEntryIdSet, Set(CallCachingEntryId(2)), singleCallCachingEntryIdSet) else Map.empty
        val ejha = createEngineJobHashingActor(
          replyTo = replyTo.ref,
          activity = activity,
          cacheLookupResponses = cacheLookupResponses)

        deathWatch watch ejha

        if (activity.readFromCache) replyTo.expectMsg(CacheMiss)
        if (activity.writeToCache) replyTo.expectMsgPF(max = 5 seconds, hint = "awaiting cache hit message") {
          case CallCacheHashes(hashes) => hashes.size should be(4)
          case x => fail(s"Cache hit anticipated! Instead got a ${x.getClass.getSimpleName}")
        }

        deathWatch.expectTerminated(ejha, 5 seconds)
      }
    }
  }
}

object EngineJobHashingActorSpec extends BackendSpec {
  import org.mockito.Mockito._

  def createEngineJobHashingActor
  (
    replyTo: ActorRef,
    activity: CallCachingActivity,
    jobDescriptor: BackendJobDescriptor = templateJobDescriptor(),
    initializationData: Option[BackendInitializationData] = None,
    fileHashingActor: Option[ActorRef] = None,
    cacheLookupResponses: Map[String, Set[CallCachingEntryId]] = Map.empty,
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
    val task = mock[Task]
    val call = mock[TaskCall]
    when(task.commandTemplateString).thenReturn("Do the stuff... now!!")
    when(task.outputs).thenReturn(List.empty)
    when(call.task).thenReturn(task)
    val workflowDescriptor = mock[BackendWorkflowDescriptor]
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, BackendJobDescriptorKey(call, None, 1), Map.empty, fqnMapToDeclarationMap(inputs))
    jobDescriptor
  }

  def standardCacheLookupResponses(commandTemplate: Set[CallCachingEntryId],
                                   inputCount: Set[CallCachingEntryId],
                                   backendName: Set[CallCachingEntryId],
                                   outputCount: Set[CallCachingEntryId]) = Map(
    "command template" -> commandTemplate,
    "input count" -> inputCount,
    "backend name" -> backendName,
    "output count" -> outputCount
  )

  def twice[A](block: Int => A) = {
    block(1)
    block(2)
  }
}
