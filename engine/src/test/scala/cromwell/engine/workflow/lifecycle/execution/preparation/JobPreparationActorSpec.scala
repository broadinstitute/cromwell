package cromwell.engine.workflow.lifecycle.execution.preparation

import akka.testkit.{ImplicitSender, TestActorRef}
import cromwell.core.callcaching.{DockerWithHash, FloatingDockerTagWithoutHash}
import cromwell.core.{LocallyQualifiedName, TestKitSuite}
import cromwell.docker.DockerHashActor.DockerHashSuccessResponse
import cromwell.docker.{DockerHashRequest, DockerHashResult, DockerImageIdentifier, DockerImageIdentifierWithoutHash}
import cromwell.engine.workflow.WorkflowDockerLookupActor.WorkflowDockerLookupFailure
import cromwell.engine.workflow.lifecycle.execution.preparation.CallPreparation.{BackendJobPreparationSucceeded, CallPreparationFailed, Start}
import cromwell.services.keyvalue.KeyValueServiceActor.{KvGet, KvKeyLookupFailed, KvPair}
import org.scalatest.{BeforeAndAfter, FlatSpecLike, Matchers}
import org.specs2.mock.Mockito
import wdl4s.Declaration
import wdl4s.values.{WdlString, WdlValue}

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.{Failure, Success}

class JobPreparationActorSpec extends TestKitSuite("JobPrepActorSpecSystem") with FlatSpecLike with Matchers with ImplicitSender with BeforeAndAfter with Mockito {

  behavior of "JobPreparationActor"

  // Generated fresh for each test case in 'before'
  var helper: JobPreparationTestHelper = null
  var inputs: Map[Declaration, WdlValue] = null

  before {
    helper = new JobPreparationTestHelper()
    inputs = Map.empty[Declaration, WdlValue]
  }

  it should "fail preparation if it can't evaluate inputs or runtime attributes" in {
    val exception = new Exception("Failed to prepare inputs/attributes - part of test flow")
    val failure = Failure(exception)
    val expectedResponse = CallPreparationFailed(helper.jobKey, exception)
    val actor = TestActorRef(helper.buildTestJobPreparationActor(null, null, null, failure, List.empty), self)
    actor ! Start
    expectMsg(expectedResponse)
    helper.workflowDockerLookupActor.expectNoMsg(100 millis)
  }

  it should "prepare successfully a job without docker attribute" in {
    val attributes = Map.empty[LocallyQualifiedName, WdlValue]
    val inputsAndAttributes = Success((inputs, attributes))
    val actor = TestActorRef(helper.buildTestJobPreparationActor(null, null, null, inputsAndAttributes, List.empty), self)
    actor ! Start
    expectMsgPF(5 seconds) {
      case success: BackendJobPreparationSucceeded =>
        success.jobDescriptor.maybeCallCachingEligible.dockerHash shouldBe None
    }
    helper.workflowDockerLookupActor.expectNoMsg(1 second)
  }

  it should "not ask for the docker hash if the attribute already contains a hash" in {
    val dockerValue = "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
    val attributes = Map(
      "docker" -> WdlString(dockerValue)
    )
    val inputsAndAttributes = Success((inputs, attributes))
    val actor = TestActorRef(helper.buildTestJobPreparationActor(null, null, null, inputsAndAttributes, List.empty), self)
    actor ! Start
    expectMsgPF(5 seconds) {
      case success: BackendJobPreparationSucceeded =>
        success.jobDescriptor.runtimeAttributes("docker").valueString shouldBe dockerValue
        success.jobDescriptor.maybeCallCachingEligible shouldBe DockerWithHash("library/ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950")
    }
    helper.workflowDockerLookupActor.expectNoMsg(1 second)
  }

  it should "lookup any requested key/value prefetches before (not) performing a docker hash lookup" in {
    val dockerValue = "ubuntu:latest"
    val attributes = Map (
      "docker" -> WdlString(dockerValue)
    )
    val hashResult = DockerHashResult("sha256", "71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950")
    val inputsAndAttributes = Success((inputs, attributes))
    val prefetchedKey1 = "hello"
    val prefetchedVal1 = KvPair(helper.scopedKeyMaker(prefetchedKey1), Some("world"))
    val prefetchedKey2 = "bonjour"
    val prefetchedVal2 = KvKeyLookupFailed(KvGet(helper.scopedKeyMaker(prefetchedKey2)))
    val prefetchedValues = Map(prefetchedKey1 -> prefetchedVal1, prefetchedKey2 -> prefetchedVal2)
    var keysToPrefetch = List(prefetchedKey1, prefetchedKey2)
    val actor = TestActorRef(helper.buildTestJobPreparationActor(1 minute, 1 minutes, List.empty, inputsAndAttributes, List(prefetchedKey1, prefetchedKey2)), self)
    actor ! Start

    def respondFromKv() = {
      helper.serviceRegistryProbe.expectMsgPF(max = 100 milliseconds) {
        case KvGet(k) if keysToPrefetch.contains(k.key) =>
          actor.tell(msg = prefetchedValues(k.key), sender = helper.serviceRegistryProbe.ref)
          keysToPrefetch = keysToPrefetch diff List(k.key)
      }
    }
    respondFromKv()
    helper.workflowDockerLookupActor.expectNoMsg(max = 100 milliseconds)
    respondFromKv()

    val req = helper.workflowDockerLookupActor.expectMsgClass(classOf[DockerHashRequest])
    helper.workflowDockerLookupActor.reply(DockerHashSuccessResponse(hashResult, req))
    expectMsgPF(5 seconds) {
      case success: BackendJobPreparationSucceeded =>
        success.jobDescriptor.prefetchedKvStoreEntries should be(Map(prefetchedKey1 -> prefetchedVal1, prefetchedKey2 -> prefetchedVal2))
    }
  }

  it should "leave the docker attribute as is and provide a DockerWithHash value" in {
    val dockerValue = "ubuntu:latest"
    val attributes = Map (
      "docker" -> WdlString(dockerValue)
    )
    val hashResult = DockerHashResult("sha256", "71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950")
    val inputsAndAttributes = Success((inputs, attributes))
    val finalValue = "library/ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
    val actor = TestActorRef(helper.buildTestJobPreparationActor(1 minute, 1 minutes, List.empty, inputsAndAttributes, List.empty), self)
    actor ! Start
    helper.workflowDockerLookupActor.expectMsgClass(classOf[DockerHashRequest])
    helper.workflowDockerLookupActor.reply(DockerHashSuccessResponse(hashResult, mock[DockerHashRequest]))
    expectMsgPF(5 seconds) {
      case success: BackendJobPreparationSucceeded =>
        success.jobDescriptor.runtimeAttributes("docker").valueString shouldBe dockerValue
        success.jobDescriptor.maybeCallCachingEligible shouldBe DockerWithHash(finalValue)
    }
  }

  it should "not provide a DockerWithHash value if it can't get the docker hash" in {
    val dockerValue = "ubuntu:latest"
    val request = DockerHashRequest(DockerImageIdentifier.fromString(dockerValue).get.asInstanceOf[DockerImageIdentifierWithoutHash])
    val attributes = Map (
      "docker" -> WdlString(dockerValue)
    )
    val inputsAndAttributes = Success((inputs, attributes))
    val actor = TestActorRef(helper.buildTestJobPreparationActor(1 minute, 1 minutes, List.empty, inputsAndAttributes, List.empty), self)
    actor ! Start
    helper.workflowDockerLookupActor.expectMsgClass(classOf[DockerHashRequest])
    helper.workflowDockerLookupActor.reply(WorkflowDockerLookupFailure(new Exception("Failed to get docker hash - part of test flow"), request))
    expectMsgPF(5 seconds) {
      case success: BackendJobPreparationSucceeded =>
        success.jobDescriptor.runtimeAttributes("docker").valueString shouldBe dockerValue
        success.jobDescriptor.maybeCallCachingEligible shouldBe FloatingDockerTagWithoutHash("library/ubuntu:latest")
    }
  }
}
