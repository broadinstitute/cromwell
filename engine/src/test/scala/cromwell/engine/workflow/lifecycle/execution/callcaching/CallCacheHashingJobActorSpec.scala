package cromwell.engine.workflow.lifecycle.execution.callcaching

import _root_.wdl.draft2.model.command.StringCommandPart
import akka.actor.{ActorRef, Props}
import akka.testkit.{TestFSMRef, TestProbe}
import cats.data.NonEmptyList
import cats.syntax.validated._
import cromwell.backend._
import cromwell.backend.standard.callcaching.StandardFileHashingActor.{FileHashResponse, SingleFileHashRequest}
import cromwell.core.TestKitSuite
import cromwell.core.callcaching.{HashingFailedMessage, _}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheHashingJobActor.{CCHJAFileHashResponse, CallCacheHashingJobActorData, CompleteFileHashingResult, HashingFiles, InitialHashingResult, NextBatchOfFileHashesRequest, NoFileHashesResult, PartialFileHashingResult, WaitingForHashFileRequest}
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CacheMiss
import cromwell.util.WomMocks
import org.scalatest.concurrent.Eventually
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpecLike, Matchers}
import wom.core.LocallyQualifiedName
import wom.graph.WomIdentifier
import wom.values.{WomInteger, WomSingleFile, WomString, WomValue}

import scala.util.control.NoStackTrace

class CallCacheHashingJobActorSpec extends TestKitSuite with FlatSpecLike with BackendSpec with Matchers with Eventually with TableDrivenPropertyChecks {
  behavior of "CallCacheReadingJobActor"

  def templateJobDescriptor(inputs: Map[LocallyQualifiedName, WomValue] = Map.empty) = {
    val task = WomMocks.mockTaskDefinition("task").copy(commandTemplateBuilder = Function.const(List(StringCommandPart("Do the stuff... now!!")).validNel))
    val call = WomMocks.mockTaskCall(WomIdentifier("call"), definition = task)
    val workflowDescriptor = mock[BackendWorkflowDescriptor]
    val runtimeAttributes = Map(
      "cpu" -> WomInteger(1),
      "memory" -> WomString("3 GB"),
      "continueOnReturnCode" -> WomInteger(0),
      "docker" -> WomString("ubuntu:latest")
    )
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, BackendJobDescriptorKey(call, None, 1), runtimeAttributes, fqnWdlMapToDeclarationMap(inputs), NoDocker, None, Map.empty)
    jobDescriptor
  }

  it should "die immediately if created without cache read actor and write to cache turned off" in {
    val parent = TestProbe()
    val testActor = TestFSMRef(new CallCacheHashingJobActor(
      templateJobDescriptor(),
      None,
      None,
      Set.empty,
      "backedName",
      Props.empty,
      DockerWithHash("ubuntu@sha256:blablablba"),
      CallCachingActivity(readWriteMode = ReadCache),
      callCachePathPrefixes = None
    ), parent.ref)
    watch(testActor)
    expectTerminated(testActor)
    parent.expectMsgClass(classOf[InitialHashingResult])
    parent.expectMsg(CacheMiss)
  }

  it should "send a correct InitialHashingResult upon starting" in {
    val parent = TestProbe()
    val inputs = Map(
      "stringInput" -> WomString("hello"),
      "fileInput" -> WomSingleFile("world")
    )
    // Do not include "memory" on purpose, even though it's in the map of runtime attributes.
    // This way we can verify that only attributes with a RuntimeAttributeDefinition are used for hashing
    // Vice versa include a "failOnStderr" definition even though it's not in the map.
    // This ensures that we still record the fact that there was no failOnStderr attribute with a "N/A" value
    val runtimeAttributeDefinitions = Set(
      RuntimeAttributeDefinition("docker", None, usedInCallCaching = true),
      RuntimeAttributeDefinition("failOnStderr", None, usedInCallCaching = true),
      RuntimeAttributeDefinition("continueOnReturnCode", Option(WomInteger(0)), usedInCallCaching = true),
      RuntimeAttributeDefinition("cpu", None, usedInCallCaching = false)
    )
    val callCacheRead = TestProbe()
    val jobDescriptor: BackendJobDescriptor = templateJobDescriptor(inputs)
    val actorUnderTest = TestFSMRef(new CallCacheHashingJobActor(
      jobDescriptor,
      Option(callCacheRead.ref),
      None,
      runtimeAttributeDefinitions,
      "backedName",
      Props.empty,
      DockerWithHash("ubuntu@sha256:blablablba"),
      CallCachingActivity(readWriteMode = ReadAndWriteCache),
      callCachePathPrefixes = None
    ), parent.ref)
    
    val expectedInitialHashes = Set(
      // md5 of Do the stuff... now
      HashResult(HashKey("command template"), HashValue("2259B15D9120F50C1BD4B2A3E2CE5A0E")),
      // md5 of backendName
      HashResult(HashKey("backend name"), HashValue("DC3D1A5AB4B8064660ADE07FFDECBFFE")),
      // md5 of 2
      HashResult(HashKey("input count"), HashValue("C81E728D9D4C2F636F067F89CC14862C")),
      // md5 of 0
      HashResult(HashKey("output count"), HashValue("CFCD208495D565EF66E7DFF9F98764DA")),
      HashResult(HashKey("runtime attribute", "failOnStderr"), HashValue("N/A")),
      // md5 of 1
      HashResult(HashKey(checkForHitOrMiss = false, "runtime attribute", "cpu"), HashValue("C4CA4238A0B923820DCC509A6F75849B")),
      // md5 of 0
      HashResult(HashKey("runtime attribute", "continueOnReturnCode"), HashValue("CFCD208495D565EF66E7DFF9F98764DA")),
      // md5 of "hello" (with quotes)
      HashResult(HashKey("input", "String stringInput"), HashValue("5DEAEE1C1332199E5B5BC7C5E4F7F0C2")),
      // md5 of ubuntu@sha256:blablablba - make sure we use the dockerWithHash and not the docker runtime attribute
      HashResult(HashKey("runtime attribute", "docker"), HashValue("C811916EA68009B0EFE0A3A86D73280E"))
    )
    val expectedAggregatedInitialHash = "F1A7118BED69B5A976A17C83FBA0D9D6"
    val expectedInitialHashResult = InitialHashingResult(expectedInitialHashes, expectedAggregatedInitialHash)
    parent.expectMsg(expectedInitialHashResult)
    callCacheRead.expectMsg(expectedInitialHashResult)
    actorUnderTest.stateName shouldBe WaitingForHashFileRequest
    actorUnderTest.stateData shouldBe CallCacheHashingJobActorData(
      List(SingleFileHashRequest(jobDescriptor.key, HashKey("input", "File fileInput"), WomSingleFile("world"), None)),
      Option(callCacheRead.ref)
    )
  }
  
  def makeCCHJA(callCacheReader: Option[ActorRef],
                testFileHashingActor: ActorRef,
                parent: ActorRef,
                writeToCache: Boolean = true,
                addFileHashMockResult: Option[(CallCacheHashingJobActorData, Option[CCHJAFileHashResponse])] = None) = {
    TestFSMRef(new CallCacheHashingJobActor(
      templateJobDescriptor(),
      callCacheReader,
      None,
      Set.empty,
      "backend",
      Props.empty,
      DockerWithHash("ubuntu@256:blablabla"),
      CallCachingActivity(readWriteMode = if (writeToCache) ReadAndWriteCache else ReadCache),
      callCachePathPrefixes = None
    ) {
      override def makeFileHashingActor() = testFileHashingActor
      override def addFileHash(hashResult: HashResult, data: CallCacheHashingJobActorData) = {
        addFileHashMockResult.getOrElse(super.addFileHash(hashResult, data))
      }
    }, parent)
  }
  
  it should "send hash file requests when receiving a NextBatchOfFileHashesRequest" in {
    val callCacheReadProbe = TestProbe()
    val fileHashingActor = TestProbe()
    
    val cchja = makeCCHJA(Option(callCacheReadProbe.ref), fileHashingActor.ref, TestProbe().ref)

    val fileHashRequest1 = SingleFileHashRequest(null, null, null, null)
    val fileHashRequest2 = SingleFileHashRequest(null, null, null, null)
    cchja.setState(
      WaitingForHashFileRequest,
      CallCacheHashingJobActorData(List(List(fileHashRequest1, fileHashRequest2)), List.empty, None)
    )

    cchja ! NextBatchOfFileHashesRequest

    fileHashingActor.expectMsg(fileHashRequest1)
    fileHashingActor.expectMsg(fileHashRequest2)
    cchja.stateName shouldBe HashingFiles
  }

  it should "send NoFileHashesResult and stop if there are no input files" in {
    val parent = TestProbe()
    val callCacheReadProbe = TestProbe()
    
    val cchja = makeCCHJA(Option(callCacheReadProbe.ref), TestProbe().ref, parent.ref)
    parent.watch(cchja)
    
    cchja.setState(
      WaitingForHashFileRequest,
      CallCacheHashingJobActorData(List.empty, List.empty, Option(callCacheReadProbe.ref))
    )

    cchja ! NextBatchOfFileHashesRequest

    callCacheReadProbe.expectMsgClass(classOf[InitialHashingResult])
    parent.expectMsgClass(classOf[InitialHashingResult])
    
    callCacheReadProbe.expectMsg(NoFileHashesResult)
    parent.expectMsg(NoFileHashesResult)
    parent.expectTerminated(cchja)
  }
  
  def selfSendNextBacthRequest(ccReader: Option[ActorRef]) = {
    val fileHashingActor = TestProbe()
    val result: PartialFileHashingResult = PartialFileHashingResult(NonEmptyList.of(mock[HashResult]))
    val fileHashRequest = SingleFileHashRequest(null, null, null, null)
    val newData = CallCacheHashingJobActorData(List(List(fileHashRequest)), List.empty, ccReader)
    // still gives a CCReader when instantiating the actor, but not in the data (above)
    // This ensures the check is done with the data and not the actor attribute, as the data will change if the ccreader dies but the actor attribute
    // will stay Some(...)
    val cchja = makeCCHJA(Option(TestProbe().ref), fileHashingActor.ref, TestProbe().ref, writeToCache = true, Option(newData -> Option(result)))
    watch(cchja)
    cchja.setState(HashingFiles)

    cchja ! FileHashResponse(mock[HashResult])

    // This proves that the ccjha keeps hashing files even though there is no ccreader requesting more hashes
    fileHashingActor.expectMsg(fileHashRequest)
    cchja.stateName shouldBe HashingFiles
  }
  
  List(
    ("send itself a NextBatchOfFileHashesRequest when a batch is complete and there is no CC Reader", None),
    ("send itself a NextBatchOfFileHashesRequest when a batch is complete and there is a CC Reader", Option(TestProbe().ref))
  ) foreach {
    case (description, ccReader) => 
      it should description in selfSendNextBacthRequest(ccReader)
  }

  it should "send FinalFileHashingResult to parent and CCReader and die" in {
    val parent = TestProbe()
    val callCacheReadProbe = TestProbe()
    List(CompleteFileHashingResult(Set(mock[HashResult]), "AggregatedFileHash"), NoFileHashesResult) foreach { result =>
      val newData = CallCacheHashingJobActorData(List.empty, List.empty, Option(callCacheReadProbe.ref))
      val cchja = makeCCHJA(Option(callCacheReadProbe.ref), TestProbe().ref, parent.ref, writeToCache = true, Option(newData -> Option(result)))

      parent.expectMsgClass(classOf[InitialHashingResult])
      callCacheReadProbe.expectMsgClass(classOf[InitialHashingResult])

      parent.watch(cchja)
      cchja.setState(HashingFiles)

      cchja ! FileHashResponse(mock[HashResult])

      callCacheReadProbe.expectMsg(result)
      parent.expectMsg(result)
      parent.expectTerminated(cchja)
    }
  }

  it should "wait for next file hash if the batch is not complete yet" in {
    val callCacheReadProbe = TestProbe()
    val parent = TestProbe()
    val newData: CallCacheHashingJobActorData = CallCacheHashingJobActorData(List.empty, List.empty, Option(callCacheReadProbe.ref))
    val cchja = makeCCHJA(Option(callCacheReadProbe.ref), TestProbe().ref, parent.ref, writeToCache = true, Option(newData -> None))

    parent.expectMsgClass(classOf[InitialHashingResult])
    callCacheReadProbe.expectMsgClass(classOf[InitialHashingResult])

    cchja.setState(HashingFiles)

    cchja ! FileHashResponse(mock[HashResult])

    callCacheReadProbe.expectNoMessage()
    parent.expectNoMessage()
    cchja.stateName shouldBe HashingFiles
  }

  it should "stop if the read actor dies and writeToCache is off" in {
    val callCacheReadProbe = TestProbe()
    val parent = TestProbe()
    val cchja = makeCCHJA(Option(callCacheReadProbe.ref), TestProbe().ref, parent.ref, writeToCache = false)
    parent.watch(cchja)
    system stop callCacheReadProbe.ref
    parent.expectMsgClass(classOf[InitialHashingResult])
    parent.expectTerminated(cchja)
  }

  it should "keep going if the read actor dies and writeToCache is on" in {
    val callCacheReadProbe = TestProbe()
    val parent = TestProbe()
    val fileHasher = TestProbe()
    val cchja = makeCCHJA(Option(callCacheReadProbe.ref), fileHasher.ref, parent.ref)
    parent.expectMsgClass(classOf[InitialHashingResult])

    val hashKey = HashKey("file")
    val fileHashRequest: SingleFileHashRequest = SingleFileHashRequest(null, hashKey, null, null)
    val data = CallCacheHashingJobActorData(List(List(fileHashRequest)), List.empty, Option(callCacheReadProbe.ref))
    
    cchja.setState(WaitingForHashFileRequest, data)
    
    system stop callCacheReadProbe.ref
    fileHasher.expectMsg(fileHashRequest)
    val result: HashResult = HashResult(hashKey, HashValue("fileHash"))
    fileHasher.reply(FileHashResponse(result))

    parent.expectMsg(CompleteFileHashingResult(Set(result), "45F27DD26834DBACBB05BBB1D651F5D1"))
  }
  
  it should "propagate HashingFailedMessage errors and die" in {
    val callCacheReadProbe = TestProbe()
    val parent = TestProbe()
    val cchja = makeCCHJA(Option(callCacheReadProbe.ref), TestProbe().ref, parent.ref)
    parent.watch(cchja)
    cchja.setState(WaitingForHashFileRequest)
    parent.expectMsgClass(classOf[InitialHashingResult])
    callCacheReadProbe.expectMsgClass(classOf[InitialHashingResult])
    
    val hashFailed = HashingFailedMessage(
      "fileName",
      new Exception("Hashing failed ! - part of test flow") with NoStackTrace
    )
    cchja ! hashFailed
    parent.expectMsg(hashFailed)
    callCacheReadProbe.expectMsg(hashFailed)
    parent.expectTerminated(cchja)
  }
  
  it should "run properly when writeToCache is ON and there is no CCRead actor" in {
    val parent = TestProbe()
    val cchja = makeCCHJA(None, TestProbe().ref, parent.ref)
    parent.watch(cchja)
    
    parent.expectMsgClass(classOf[InitialHashingResult])
    parent.expectMsg(NoFileHashesResult)
    parent.expectTerminated(cchja)
  }
}
