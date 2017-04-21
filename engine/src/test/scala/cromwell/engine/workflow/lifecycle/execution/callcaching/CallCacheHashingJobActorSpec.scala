package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.actor.{ActorRef, Props}
import akka.testkit.{TestFSMRef, TestProbe}
import cromwell.backend._
import cromwell.backend.standard.callcaching.StandardFileHashingActor.{FileHashResponse, SingleFileHashRequest}
import cromwell.core.callcaching.{HashingFailedMessage, _}
import cromwell.core.{LocallyQualifiedName, TestKitSuite}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheHashingJobActor.{CCHJAFileHashResponse, CallCacheHashingJobActorData, CompleteFileHashingResult, HashingFiles, InitialHashingResult, NextBatchOfFileHashesRequest, NoFileHashesResult, PartialFileHashingResult, WaitingForHashFileRequest}
import org.mockito.Mockito._
import org.scalatest.concurrent.Eventually
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpecLike, Matchers}
import wdl4s.values.{WdlFile, WdlInteger, WdlString, WdlValue}
import wdl4s.{Task, TaskCall}

class CallCacheHashingJobActorSpec extends TestKitSuite with FlatSpecLike with BackendSpec with Matchers with Eventually with TableDrivenPropertyChecks {
  behavior of "CallCacheReadingJobActor"

  def templateJobDescriptor(inputs: Map[LocallyQualifiedName, WdlValue] = Map.empty) = {
    val task = mock[Task]
    val call = mock[TaskCall]
    when(task.commandTemplateString).thenReturn("Do the stuff... now!!")
    when(task.outputs).thenReturn(List.empty)
    when(call.task).thenReturn(task)
    val workflowDescriptor = mock[BackendWorkflowDescriptor]
    val runtimeAttributes = Map(
      "cpu" -> WdlInteger(1),
      "memory" -> WdlString("3 GB"),
      "continueOnReturnCode" -> WdlInteger(0)
    )
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, BackendJobDescriptorKey(call, None, 1), runtimeAttributes, fqnMapToDeclarationMap(inputs), CallCachingEligible, Map.empty)
    jobDescriptor
  }

  it should "die immediately if created without cache read actor and write to cache turned off" in {
    val testActor = TestFSMRef(new CallCacheHashingJobActor(
      templateJobDescriptor(),
      None,
      None,
      Set.empty,
      "backedName",
      Props.empty,
      false
    ))
    watch(testActor)
    expectTerminated(testActor)
  }

  it should "send a correct InitialHashingResult upon starting" in {
    val parent = TestProbe()
    val inputs = Map(
      "stringInput" -> WdlString("hello"),
      "fileInput" -> WdlFile("world")
    )
    // Do not include "memory" on purpose, even though it's in the map of runtime attributes.
    // This way we can verify that only attributes with a RuntimeAttributeDefinition are used for hashing
    // Vice versa include a "docker" definition even though it's not in the map.
    // This ensures that we still record the fact that there was no docker attribute with a "N/A" value
    val runtimeAttributeDefinitions = Set(
      RuntimeAttributeDefinition("docker", None, usedInCallCaching = true),
      RuntimeAttributeDefinition("continueOnReturnCode", Option(WdlInteger(0)), usedInCallCaching = true),
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
      false
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
      HashResult(HashKey("runtime attribute: docker"), HashValue("N/A")),
      // md5 of 1
      HashResult(HashKey("runtime attribute: cpu", checkForHitOrMiss = false), HashValue("C4CA4238A0B923820DCC509A6F75849B")),
      // md5 of 0
      HashResult(HashKey("runtime attribute: continueOnReturnCode"), HashValue("CFCD208495D565EF66E7DFF9F98764DA")),
      // md5 of "hello" (with quotes)
      HashResult(HashKey("input: String stringInput"), HashValue("5DEAEE1C1332199E5B5BC7C5E4F7F0C2"))
    )
    val expectedAggregatedInitialHash = "67A8FA54B6D3A57D40AA6FE39C6C3992"
    val expectedInitialHashResult = InitialHashingResult(expectedInitialHashes, expectedAggregatedInitialHash)
    parent.expectMsg(expectedInitialHashResult)
    callCacheRead.expectMsg(expectedInitialHashResult)
    actorUnderTest.stateName shouldBe WaitingForHashFileRequest
    actorUnderTest.stateData shouldBe CallCacheHashingJobActorData(
      List(SingleFileHashRequest(jobDescriptor.key, HashKey("input: File fileInput"), WdlFile("world"), None)),
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
      writeToCache = writeToCache
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
    // self acts as the parent
    val cchja = makeCCHJA(Option(callCacheReadProbe.ref), fileHashingActor.ref, TestProbe().ref)

    val fileHashRequest1 = mock[SingleFileHashRequest]
    val fileHashRequest2 = mock[SingleFileHashRequest]
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
    // self acts as the parent
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

  it should "send the PartialFileHashingResult when a batch is complete" in {
    val callCacheReadProbe = TestProbe()
    val hashResults = Set(mock[HashResult])
    
    val result: PartialFileHashingResult = PartialFileHashingResult(hashResults)
    val newData: CallCacheHashingJobActorData = CallCacheHashingJobActorData(List.empty, List.empty, Option(callCacheReadProbe.ref))
    // self acts as the parent
    val cchja = makeCCHJA(Option(callCacheReadProbe.ref), TestProbe().ref, TestProbe().ref, writeToCache = true, Option(newData -> Option(result)))

    cchja.setState(HashingFiles)

    cchja ! FileHashResponse(mock[HashResult])

    callCacheReadProbe.expectMsgClass(classOf[InitialHashingResult])
    callCacheReadProbe.expectMsg(result)
    cchja.stateName shouldBe WaitingForHashFileRequest
    cchja.stateData shouldBe newData
  }

  it should "send itself a NextBatchOfFileHashesRequest when a batch is complete and there is no CCReader" in {
    val fileHashingActor = TestProbe()
    val result: PartialFileHashingResult = PartialFileHashingResult(Set(mock[HashResult]))
    val fileHashRequest = mock[SingleFileHashRequest]
    val newData = CallCacheHashingJobActorData(List(List(fileHashRequest)), List.empty, None)
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

    callCacheReadProbe.expectNoMsg()
    parent.expectNoMsg()
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
    val cchja = makeCCHJA(Option(callCacheReadProbe.ref), TestProbe().ref, parent.ref)
    parent.expectMsgClass(classOf[InitialHashingResult])
    
    cchja.setState(WaitingForHashFileRequest)
    system stop callCacheReadProbe.ref

    parent.expectMsg(NoFileHashesResult)
  }
  
  it should "propagate HashingFailedMessage errors and die" in {
    val callCacheReadProbe = TestProbe()
    val parent = TestProbe()
    val cchja = makeCCHJA(Option(callCacheReadProbe.ref), TestProbe().ref, parent.ref)
    parent.watch(cchja)
    cchja.setState(WaitingForHashFileRequest)
    parent.expectMsgClass(classOf[InitialHashingResult])
    callCacheReadProbe.expectMsgClass(classOf[InitialHashingResult])
    
    val hashFailed = HashingFailedMessage("fileName", new Exception("Hashing failed ! - part of test flow"))
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
