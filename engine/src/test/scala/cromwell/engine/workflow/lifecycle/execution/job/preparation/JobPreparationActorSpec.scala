package cromwell.engine.workflow.lifecycle.execution.job.preparation

import akka.testkit.{ImplicitSender, TestActorRef}
import cats.syntax.validated._
import common.mock.MockSugar
import cromwell.backend.BackendJobDescriptorKey
import cromwell.core.TestKitSuite
import cromwell.core.callcaching.{DockerWithHash, FloatingDockerTagWithoutHash}
import cromwell.docker.DockerInfoActor.{DockerInformation, DockerInfoSuccessResponse, DockerSize}
import cromwell.docker._
import cromwell.engine.workflow.WorkflowDockerLookupActor.WorkflowDockerLookupFailure
import cromwell.engine.workflow.lifecycle.execution.job.preparation.CallPreparation.{
  BackendJobPreparationSucceeded,
  CallPreparationFailed,
  Start
}
import cromwell.engine.workflow.lifecycle.execution.stores.ValueStore
import cromwell.services.keyvalue.KeyValueServiceActor.{KvGet, KvKeyLookupFailed, KvPair}
import org.scalatest.BeforeAndAfter
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import wom.{RuntimeAttributes, RuntimeAttributesKeys}
import wom.callable.Callable.InputDefinition
import wom.callable.CommandTaskDefinition
import wom.core.LocallyQualifiedName
import wom.expression.ValueAsAnExpression
import wom.graph.{CommandCallNode, WomIdentifier}
import wom.types.{WomArrayType, WomStringType}
import wom.values.{WomArray, WomString, WomValue}

import scala.concurrent.duration._
import scala.language.postfixOps
import scala.util.control.NoStackTrace

class JobPreparationActorSpec
    extends TestKitSuite
    with AnyFlatSpecLike
    with Matchers
    with ImplicitSender
    with BeforeAndAfter
    with MockSugar {

  behavior of "JobPreparationActor"

  // Generated fresh for each test case in 'before'
  var helper: JobPreparationTestHelper = _
  var inputs: Map[InputDefinition, WomValue] = _

  before {
    helper = new JobPreparationTestHelper()
    inputs = Map.empty[InputDefinition, WomValue]
  }

  it should "fail preparation if it can't evaluate inputs or runtime attributes" in {
    val error = "Failed to prepare inputs/attributes - part of test flow"
    val actor =
      TestActorRef(helper.buildTestJobPreparationActor(null, null, null, Option(error.invalidNel), List.empty), self)
    actor ! Start(ValueStore.empty)
    expectMsgPF(1.second) { case CallPreparationFailed(_, ex) =>
      ex.getMessage shouldBe "Call input and runtime attributes evaluation failed for JobPreparationSpec_call:\nFailed to prepare inputs/attributes - part of test flow"
    }
    helper.workflowDockerLookupActor.expectNoMessage(100 millis)
  }

  it should "prepare successfully a job without docker attribute" in {
    val attributes = Map.empty[LocallyQualifiedName, WomValue]
    val inputsAndAttributes = (inputs, attributes).validNel
    val actor =
      TestActorRef(helper.buildTestJobPreparationActor(null, null, null, Option(inputsAndAttributes), List.empty), self)
    actor ! Start(ValueStore.empty)
    expectMsgPF(5 seconds) { case success: BackendJobPreparationSucceeded =>
      success.jobDescriptor.maybeCallCachingEligible.dockerHash shouldBe None
    }
    helper.workflowDockerLookupActor.expectNoMessage(1 second)
  }

  it should "ask for the docker info even if the attribute already contains a hash" in {
    val dockerValue = "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
    val attributes = Map(
      "docker" -> WomString(dockerValue)
    )
    val inputsAndAttributes = (inputs, attributes).validNel
    val actor =
      TestActorRef(helper.buildTestJobPreparationActor(null, null, null, Option(inputsAndAttributes), List.empty), self)
    actor ! Start(ValueStore.empty)
    helper.workflowDockerLookupActor.expectMsgClass(classOf[DockerInfoRequest])
    actor ! DockerInfoSuccessResponse(
      DockerInformation(
        DockerHashResult("sha256", "71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"),
        Option(DockerSize(100))
      ),
      null
    )
    expectMsgPF(5 seconds) { case success: BackendJobPreparationSucceeded =>
      success.jobDescriptor.runtimeAttributes("docker").valueString shouldBe dockerValue
      success.jobDescriptor.maybeCallCachingEligible shouldBe DockerWithHash(
        "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
      )
      success.jobDescriptor.dockerSize shouldBe Option(DockerSize(100))
    }
  }

  it should "lookup any requested key/value prefetches after (not) performing a docker hash lookup" in {
    val dockerValue = "ubuntu:latest"
    val attributes = Map(
      "docker" -> WomString(dockerValue)
    )
    val hashResult = DockerHashResult("sha256", "71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950")
    val inputsAndAttributes = (inputs, attributes).validNel
    val prefetchedKey1 = "hello"
    val prefetchedVal1 = KvPair(helper.scopedKeyMaker(prefetchedKey1), "world")
    val prefetchedKey2 = "bonjour"
    val prefetchedVal2 = KvKeyLookupFailed(KvGet(helper.scopedKeyMaker(prefetchedKey2)))
    val prefetchedValues = Map(prefetchedKey1 -> prefetchedVal1, prefetchedKey2 -> prefetchedVal2)
    var keysToPrefetch = List(prefetchedKey1, prefetchedKey2)
    val actor = TestActorRef(helper.buildTestJobPreparationActor(1 minute,
                                                                 1 minutes,
                                                                 List.empty,
                                                                 Option(inputsAndAttributes),
                                                                 List(prefetchedKey1, prefetchedKey2)
                             ),
                             self
    )
    actor ! Start(ValueStore.empty)

    val req = helper.workflowDockerLookupActor.expectMsgClass(classOf[DockerInfoRequest])
    helper.workflowDockerLookupActor.reply(DockerInfoSuccessResponse(DockerInformation(hashResult, None), req))

    def respondFromKv(): Unit =
      helper.serviceRegistryProbe.expectMsgPF(max = 100 milliseconds) {
        case KvGet(k) if keysToPrefetch.contains(k.key) =>
          actor.tell(msg = prefetchedValues(k.key), sender = helper.serviceRegistryProbe.ref)
          keysToPrefetch = keysToPrefetch diff List(k.key)
      }
    respondFromKv()
    helper.workflowDockerLookupActor.expectNoMessage(max = 100 milliseconds)
    respondFromKv()

    expectMsgPF(5 seconds) { case success: BackendJobPreparationSucceeded =>
      success.jobDescriptor.prefetchedKvStoreEntries should be(
        Map(prefetchedKey1 -> prefetchedVal1, prefetchedKey2 -> prefetchedVal2)
      )
    }
  }

  it should "leave the docker attribute as is and provide a DockerWithHash value" in {
    val dockerValue = "ubuntu:latest"
    val attributes = Map(
      "docker" -> WomString(dockerValue)
    )
    val hashResult = DockerHashResult("sha256", "71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950")
    val inputsAndAttributes = (inputs, attributes).validNel
    val finalValue = "ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
    val actor = TestActorRef(
      helper.buildTestJobPreparationActor(1 minute, 1 minutes, List.empty, Option(inputsAndAttributes), List.empty),
      self
    )
    actor ! Start(ValueStore.empty)
    helper.workflowDockerLookupActor.expectMsgClass(classOf[DockerInfoRequest])
    helper.workflowDockerLookupActor.reply(
      DockerInfoSuccessResponse(DockerInformation(hashResult, None), mock[DockerInfoRequest])
    )
    expectMsgPF(5 seconds) { case success: BackendJobPreparationSucceeded =>
      success.jobDescriptor.runtimeAttributes("docker").valueString shouldBe dockerValue
      success.jobDescriptor.maybeCallCachingEligible shouldBe DockerWithHash(finalValue)
    }
  }

  it should "not provide a DockerWithHash value if it can't get the docker hash" in {
    val dockerValue = "ubuntu:latest"
    val request = DockerInfoRequest(
      DockerImageIdentifier.fromString(dockerValue).get.asInstanceOf[DockerImageIdentifierWithoutHash]
    )
    val attributes = Map(
      "docker" -> WomString(dockerValue)
    )
    val inputsAndAttributes = (inputs, attributes).validNel
    val actor = TestActorRef(
      helper.buildTestJobPreparationActor(1 minute, 1 minutes, List.empty, Option(inputsAndAttributes), List.empty),
      self
    )
    actor ! Start(ValueStore.empty)
    helper.workflowDockerLookupActor.expectMsgClass(classOf[DockerInfoRequest])
    helper.workflowDockerLookupActor.reply(
      WorkflowDockerLookupFailure(
        new Exception("Failed to get docker hash - part of test flow") with NoStackTrace,
        request
      )
    )
    expectMsgPF(5 seconds) { case success: BackendJobPreparationSucceeded =>
      success.jobDescriptor.runtimeAttributes("docker").valueString shouldBe dockerValue
      success.jobDescriptor.maybeCallCachingEligible shouldBe FloatingDockerTagWithoutHash("ubuntu:latest")
    }
  }

  it should "correctly apply DockerMirroring to a Dockerhub image defined in docker attr" in {
    // Create a mock job with the desired runtime attribute
    val callable: CommandTaskDefinition = mock[CommandTaskDefinition]
    callable.runtimeAttributes returns RuntimeAttributes(
      Map("docker" -> ValueAsAnExpression(WomString("ubuntu:latest")))
    )
    val call: CommandCallNode =
      CommandCallNode(WomIdentifier("JobPreparationSpec_call"), callable, Set.empty, List.empty, Set.empty, null, None)
    val mockJobKey: BackendJobDescriptorKey = BackendJobDescriptorKey(call, None, 1)

    val dockerMirroring = DockerMirroring(List(DockerHubMirror("my.mirror.io")))
    val hashResult = DockerHashResult("sha256", "71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950")
    val mirroredValue = "my.mirror.io/library/ubuntu:latest"
    val finalValue =
      "my.mirror.io/library/ubuntu@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
    val actor = TestActorRef(
      helper.buildTestJobPreparationActor(1 minute,
                                          1 minutes,
                                          List.empty,
                                          None,
                                          List.empty,
                                          mockJobKey,
                                          Option(dockerMirroring)
      ),
      self
    )
    actor ! Start(ValueStore.empty)
    helper.workflowDockerLookupActor.expectMsgPF(5 seconds) { case success: DockerInfoRequest =>
      success.dockerImageID.fullName shouldBe mirroredValue
    }
    helper.workflowDockerLookupActor.reply(
      DockerInfoSuccessResponse(DockerInformation(hashResult, None), mock[DockerInfoRequest])
    )
    expectMsgPF(5 seconds) { case success: BackendJobPreparationSucceeded =>
      success.jobDescriptor.runtimeAttributes("docker").valueString shouldBe mirroredValue
      success.jobDescriptor.maybeCallCachingEligible shouldBe DockerWithHash(finalValue)
    }
  }

  it should "correctly apply DockerMirroring to a Dockerhub image defined in container attr" in {
    // Create a mock job with the desired runtime attribute
    val callable: CommandTaskDefinition = mock[CommandTaskDefinition]
    callable.runtimeAttributes returns RuntimeAttributes(
      Map(
        "container" -> ValueAsAnExpression(
          WomArray(WomArrayType(WomStringType), Seq(WomString("alpine:latest"), WomString("ubuntu:latest")))
        )
      )
    )
    val call: CommandCallNode =
      CommandCallNode(WomIdentifier("JobPreparationSpec_call"), callable, Set.empty, List.empty, Set.empty, null, None)
    val mockJobKey: BackendJobDescriptorKey = BackendJobDescriptorKey(call, None, 1)

    val dockerMirroring = DockerMirroring(List(DockerHubMirror("my.mirror.io")))
    val hashResult = DockerHashResult("sha256", "71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950")
    val mirroredUbuntuValue = "my.mirror.io/library/ubuntu:latest"
    val mirroredAlpineValue = "my.mirror.io/library/alpine:latest"
    val finalValue =
      "my.mirror.io/library/alpine@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
    val actor = TestActorRef(
      helper.buildTestJobPreparationActor(1 minute,
                                          1 minutes,
                                          List.empty,
                                          None,
                                          List.empty,
                                          mockJobKey,
                                          Option(dockerMirroring)
      ),
      self
    )
    actor ! Start(ValueStore.empty)
    helper.workflowDockerLookupActor.expectMsgPF(5 seconds) { case success: DockerInfoRequest =>
      success.dockerImageID.fullName shouldBe mirroredAlpineValue
    }
    helper.workflowDockerLookupActor.reply(
      DockerInfoSuccessResponse(DockerInformation(hashResult, None), mock[DockerInfoRequest])
    )
    expectMsgPF(5 seconds) { case success: BackendJobPreparationSucceeded =>
      success.jobDescriptor
        .runtimeAttributes("container")
        .toString shouldBe s"[\"${mirroredAlpineValue}\", \"${mirroredUbuntuValue}\"]"
      success.jobDescriptor.maybeCallCachingEligible shouldBe DockerWithHash(finalValue)
    }
  }

  it should "not apply DockerMirroring to a non-DockerHub image" in {
    val dockerValue = "gcr.io/broad-dsde-cromwell-dev/cromwell-drs-localizer:latest"

    // Create a mock job with the desired runtime attribute
    val callable: CommandTaskDefinition = mock[CommandTaskDefinition]
    callable.runtimeAttributes returns RuntimeAttributes(
      Map("docker" -> ValueAsAnExpression(WomString(dockerValue)))
    )
    val call: CommandCallNode =
      CommandCallNode(WomIdentifier("JobPreparationSpec_call"), callable, Set.empty, List.empty, Set.empty, null, None)
    val mockJobKey: BackendJobDescriptorKey = BackendJobDescriptorKey(call, None, 1)

    val dockerMirroring = DockerMirroring(List(DockerHubMirror("my.mirror.io")))
    val hashResult = DockerHashResult("sha256", "71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950")
    val finalValue =
      "gcr.io/broad-dsde-cromwell-dev/cromwell-drs-localizer@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
    val actor = TestActorRef(
      helper.buildTestJobPreparationActor(1 minute,
                                          1 minutes,
                                          List.empty,
                                          None,
                                          List.empty,
                                          mockJobKey,
                                          Option(dockerMirroring)
      ),
      self
    )
    actor ! Start(ValueStore.empty)
    helper.workflowDockerLookupActor.expectMsgClass(classOf[DockerInfoRequest])
    helper.workflowDockerLookupActor.reply(
      DockerInfoSuccessResponse(DockerInformation(hashResult, None), mock[DockerInfoRequest])
    )
    expectMsgPF(5 seconds) { case success: BackendJobPreparationSucceeded =>
      success.jobDescriptor.runtimeAttributes("docker").valueString shouldBe dockerValue
      success.jobDescriptor.maybeCallCachingEligible shouldBe DockerWithHash(finalValue)
    }
  }

  it should "prefer `container` over `docker` as image source when both are present" in {
    // Create a mock job with the desired runtime attribute
    val callable: CommandTaskDefinition = mock[CommandTaskDefinition]
    callable.runtimeAttributes returns RuntimeAttributes(
      Map(
        "docker" -> ValueAsAnExpression(WomString("ubuntu:latest")),
        "container" -> ValueAsAnExpression(WomArray(WomArrayType(WomStringType), Seq(WomString("alpine:latest"))))
      )
    )
    val call: CommandCallNode =
      CommandCallNode(WomIdentifier("JobPreparationSpec_call"), callable, Set.empty, List.empty, Set.empty, null, None)
    val mockJobKey: BackendJobDescriptorKey = BackendJobDescriptorKey(call, None, 1)

    val hashResult = DockerHashResult("sha256", "71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950")
    val finalValue =
      "alpine@sha256:71cd81252a3563a03ad8daee81047b62ab5d892ebbfbf71cf53415f29c130950"
    val actor = TestActorRef(
      helper.buildTestJobPreparationActor(1 minute, 1 minutes, List.empty, None, List.empty, mockJobKey, None),
      self
    )
    actor ! Start(ValueStore.empty)
    helper.workflowDockerLookupActor.expectMsgPF(5 seconds) { case success: DockerInfoRequest =>
      success.dockerImageID.fullName shouldBe "alpine:latest"
    }
    helper.workflowDockerLookupActor.reply(
      DockerInfoSuccessResponse(DockerInformation(hashResult, None), mock[DockerInfoRequest])
    )
    expectMsgPF(5 seconds) { case success: BackendJobPreparationSucceeded =>
      success.jobDescriptor.maybeCallCachingEligible shouldBe DockerWithHash(finalValue)
    }
  }

  it should "lookup MemoryMultiplier key/value if available and accordingly update runtime attributes" in {
    val attributes = Map(
      "memory" -> WomString("1.1 GB")
    )
    val inputsAndAttributes = (inputs, attributes).validNel
    val prefetchedKey = "MemoryMultiplier"
    val prefetchedVal = KvPair(helper.scopedKeyMaker(prefetchedKey), "1.1")
    val prefetchedValues = Map(prefetchedKey -> prefetchedVal)
    var keysToPrefetch = List(prefetchedKey)
    val actor = TestActorRef(
      helper.buildTestJobPreparationActor(1 minute,
                                          1 minutes,
                                          List.empty,
                                          Option(inputsAndAttributes),
                                          List(prefetchedKey)
      ),
      self
    )
    actor ! Start(ValueStore.empty)

    helper.serviceRegistryProbe.expectMsgPF(max = 100 milliseconds) {
      case KvGet(k) if keysToPrefetch.contains(k.key) =>
        actor.tell(msg = prefetchedValues(k.key), sender = helper.serviceRegistryProbe.ref)
        keysToPrefetch = keysToPrefetch diff List(k.key)
    }

    expectMsgPF(5 seconds) { case success: BackendJobPreparationSucceeded =>
      success.jobDescriptor.prefetchedKvStoreEntries should be(Map(prefetchedKey -> prefetchedVal))
      success.jobDescriptor.runtimeAttributes(RuntimeAttributesKeys.MemoryKey) shouldBe WomString(
        "1.2100000000000002 GB"
      )
    }
  }

  it should "mimic 100 memory retries and lookup MemoryMultiplier key/value to update runtime attributes" in {
    val prefetchedKey = "MemoryMultiplier"
    val retryFactor = 1.1
    val taskMemory = 1.0
    val attributes = Map("memory" -> WomString(s"$taskMemory GB"))
    val inputsAndAttributes = (inputs, attributes).validNel

    var previousMultiplier = 1.0

    // mimic failure and retry task 100 times with more memory. At the end of this loop, the task memory would
    // have been 1.1 ^ 100 = 13780.612339822379 GB. Note: This is just a hypothetical situation. In reality, PAPI would
    // fail the task because a VM with such a large amount of memory might not exist
    for (_ <- 1 to 100) {
      // since the task has failed, increase it's multiplier
      val nextMultiplier = previousMultiplier * retryFactor
      previousMultiplier = nextMultiplier

      val prefetchedVal = KvPair(helper.scopedKeyMaker(prefetchedKey), nextMultiplier.toString)
      val prefetchedValues = Map(prefetchedKey -> prefetchedVal)
      var keysToPrefetch = List(prefetchedKey)

      val actor = TestActorRef(
        helper.buildTestJobPreparationActor(1 minute,
                                            1 minutes,
                                            List.empty,
                                            Option(inputsAndAttributes),
                                            List(prefetchedKey)
        ),
        self
      )
      actor ! Start(ValueStore.empty)

      helper.serviceRegistryProbe.expectMsgPF(max = 100 milliseconds) {
        case KvGet(k) if keysToPrefetch.contains(k.key) =>
          actor.tell(msg = prefetchedValues(k.key), sender = helper.serviceRegistryProbe.ref)
          keysToPrefetch = keysToPrefetch diff List(k.key)
      }

      expectMsgPF(5 seconds) { case success: BackendJobPreparationSucceeded =>
        success.jobDescriptor.prefetchedKvStoreEntries should be(Map(prefetchedKey -> prefetchedVal))
        success.jobDescriptor.runtimeAttributes(RuntimeAttributesKeys.MemoryKey) shouldBe WomString(
          s"${taskMemory * nextMultiplier} GB"
        )
      }
    }
  }
}
