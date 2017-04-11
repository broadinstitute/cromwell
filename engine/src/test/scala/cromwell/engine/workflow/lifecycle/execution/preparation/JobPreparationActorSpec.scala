package cromwell.engine.workflow.lifecycle.execution.preparation

import akka.testkit.{ImplicitSender, TestActorRef}
import cromwell.core.actor.StreamIntegration.BackPressure
import cromwell.core.callcaching.docker.DockerHashActor.{DockerHashFailedResponse, DockerHashResponseSuccess}
import cromwell.core.callcaching.docker.{DockerHashRequest, DockerHashResult, DockerImageIdentifier, DockerImageIdentifierWithoutHash}
import cromwell.core.callcaching.{CallCachingEligible, CallCachingIneligible}
import cromwell.core.{LocallyQualifiedName, TestKitSuite}
import cromwell.engine.workflow.lifecycle.execution.preparation.CallPreparation.{BackendJobPreparationSucceeded, CallPreparationFailed, Start}
import cromwell.services.keyvalue.KeyValueServiceActor.{KvGet, KvKeyLookupFailed, KvPair}
import org.scalatest.BeforeAndAfter
import org.scalatest.{FlatSpecLike, Matchers}
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
    helper.dockerHashingActor.expectNoMsg(100 millis)
  }

  it should "prepare successfully a job without docker attribute" in {
    val attributes = Map.empty[LocallyQualifiedName, WdlValue]
    val inputsAndAttributes = Success((inputs, attributes))
    val actor = TestActorRef(helper.buildTestJobPreparationActor(null, null, null, inputsAndAttributes, List.empty), self)
    actor ! Start
    expectMsgPF(5 seconds) {
      case success: BackendJobPreparationSucceeded =>
        success.jobDescriptor.callCachingEligibility shouldBe CallCachingEligible
    }
    helper.dockerHashingActor.expectNoMsg(1 second)
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
        success.jobDescriptor.callCachingEligibility shouldBe CallCachingEligible
    }
    helper.dockerHashingActor.expectNoMsg(1 second)
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
    helper.dockerHashingActor.expectNoMsg(max = 100 milliseconds)
    respondFromKv()

    val req = helper.dockerHashingActor.expectMsgClass(classOf[DockerHashRequest])
    helper.dockerHashingActor.reply(DockerHashResponseSuccess(hashResult, req))
    expectMsgPF(5 seconds) {
      case success: BackendJobPreparationSucceeded =>
        success.jobDescriptor.prefetchedKvStoreEntries should be(Map(prefetchedKey1 -> prefetchedVal1, prefetchedKey2 -> prefetchedVal2))
    }
  }

  it should "replace the tag with a docker hash if the attribute has a floating tag, and make it ineligible for caching" in {
    val dockerValue = "ubuntu:latest"
    val attributes = Map (
      "docker" -> WdlString(dockerValue)
    )
    val hashResult = DockerHashResult("sha256", "71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950")
    val inputsAndAttributes = Success((inputs, attributes))
    val finalValue = "library/ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
    val actor = TestActorRef(helper.buildTestJobPreparationActor(1 minute, 1 minutes, List.empty, inputsAndAttributes, List.empty), self)
    actor ! Start
    helper.dockerHashingActor.expectMsgClass(classOf[DockerHashRequest])
    helper.dockerHashingActor.reply(DockerHashResponseSuccess(hashResult, mock[DockerHashRequest]))
    expectMsgPF(5 seconds) {
      case success: BackendJobPreparationSucceeded =>
        success.jobDescriptor.runtimeAttributes("docker").valueString shouldBe finalValue
        success.jobDescriptor.callCachingEligibility.isInstanceOf[CallCachingIneligible] shouldBe true
        success.jobDescriptor.callCachingEligibility.asInstanceOf[CallCachingIneligible].message shouldBe
          s"""You are using a floating docker tag in this task. Cromwell does not consider tasks with floating tags to be eligible for call caching.
              |If you want this task to be eligible for call caching in the future, use a docker runtime attribute with a digest instead.
              |This is the exact docker image that was used for this job: $finalValue
              |You can replace the docker runtime attribute in your task with the above value to make this task eligible for call caching.""".stripMargin
    }
  }

  it should "make a job ineligible for caching if it can't get the docker hash" in {
    val dockerValue = "ubuntu:latest"
    val request = DockerHashRequest(DockerImageIdentifier.fromString(dockerValue).get.asInstanceOf[DockerImageIdentifierWithoutHash])
    val attributes = Map (
      "docker" -> WdlString(dockerValue)
    )
    val inputsAndAttributes = Success((inputs, attributes))
    val actor = TestActorRef(helper.buildTestJobPreparationActor(1 minute, 1 minutes, List.empty, inputsAndAttributes, List.empty), self)
    actor ! Start
    helper.dockerHashingActor.expectMsgClass(classOf[DockerHashRequest])
    helper.dockerHashingActor.reply(DockerHashFailedResponse(new Exception("Failed to get docker hash - part of test flow"), request))
    expectMsgPF(5 seconds) {
      case success: BackendJobPreparationSucceeded =>
        success.jobDescriptor.runtimeAttributes("docker").valueString shouldBe dockerValue
        success.jobDescriptor.callCachingEligibility.isInstanceOf[CallCachingIneligible] shouldBe true
        success.jobDescriptor.callCachingEligibility.asInstanceOf[CallCachingIneligible].message shouldBe
          s"""You are using a floating docker tag in this task. Cromwell does not consider tasks with floating tags to be eligible for call caching.
              |If you want this task to be eligible for call caching in the future, use a docker runtime attribute with a digest instead.
              |Cromwell attempted to retrieve the current hash for this docker image but failed.
              |This is not necessarily a cause for concern as Cromwell is currently only able to retrieve hashes for Dockerhub and GCR images.
              |The job will be dispatched to the appropriate backend that will attempt to run it.""".stripMargin
    }
  }
  
  it should "wait and resubmit the docker request when it gets a backpressure message" in {
    val dockerValue = "ubuntu:latest"
    val dockerId = DockerImageIdentifier.fromString(dockerValue).get.asInstanceOf[DockerImageIdentifierWithoutHash]
    val attributes = Map (
      "docker" -> WdlString(dockerValue)
    )
    val inputsAndAttributes = Success((inputs, attributes))
    val backpresureWaitTime = 2 seconds
    val actor = TestActorRef(helper.buildTestJobPreparationActor(backpresureWaitTime, 1 minute, List.empty, inputsAndAttributes, List.empty), self)
    val request = DockerHashRequest(dockerId)
    actor ! Start
    helper.dockerHashingActor.expectMsg(request)
    helper.dockerHashingActor.reply(BackPressure(request))
    // Give a couple of seconds of margin to account for test latency etc...
    helper.dockerHashingActor.expectMsg(backpresureWaitTime.+(2 seconds), request)
  }

  it should "time out if no answer is received from the docker hash actor" in {
    val dockerValue = "ubuntu:latest"
    val dockerId = DockerImageIdentifier.fromString(dockerValue).get.asInstanceOf[DockerImageIdentifierWithoutHash]
    val attributes = Map (
      "docker" -> WdlString(dockerValue)
    )
    val inputsAndAttributes = Success((inputs, attributes))
    val requestTimeout = 1 second
    val actor = TestActorRef(helper.buildTestJobPreparationActor(1 minute, requestTimeout, List.empty, inputsAndAttributes, List.empty), self)
    val request = DockerHashRequest(dockerId)
    actor ! Start
    
    within(10 seconds) {
      helper.dockerHashingActor.expectMsg(request)
    }

    expectMsgPF(5 seconds) {
      case success: BackendJobPreparationSucceeded =>
        success.jobDescriptor.runtimeAttributes("docker").valueString shouldBe dockerValue
        success.jobDescriptor.callCachingEligibility.isInstanceOf[CallCachingIneligible] shouldBe true
        success.jobDescriptor.callCachingEligibility.asInstanceOf[CallCachingIneligible].message shouldBe
          s"""You are using a floating docker tag in this task. Cromwell does not consider tasks with floating tags to be eligible for call caching.
              |If you want this task to be eligible for call caching in the future, use a docker runtime attribute with a digest instead.
              |Cromwell attempted to retrieve the current hash for this docker image but failed.
              |This is not necessarily a cause for concern as Cromwell is currently only able to retrieve hashes for Dockerhub and GCR images.
              |The job will be dispatched to the appropriate backend that will attempt to run it.""".stripMargin
    }
  }
  
}
