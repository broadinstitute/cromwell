package cromwell.backend.google.pipelines.common

import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.refineV
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

import scala.concurrent.duration._
import scala.language.postfixOps

class PipelinesApiBackendLifecycleActorFactorySpec extends AnyFlatSpecLike with Matchers with TableDrivenPropertyChecks {

  "PipelinesApiBackendLifecycleActorFactory" should "robustly build configuration attributes" in {

    val attributes = new PipelinesApiConfigurationAttributes(
      project = "project",
      computeServiceAccount = "computeServiceAccount",
      auths = null,
      restrictMetadataAccess = true,
      enableFuse = true,
      executionBucket = "executionBucket",
      endpointUrl = null,
      location = "location",
      maxPollingInterval = 0,
      qps = refineV[Positive](1).right.get,
      cacheHitDuplicationStrategy = null,
      requestWorkers = refineV[Positive](1).right.get,
      pipelineTimeout = 1 second,
      logFlushPeriod = Option(1 second),
      gcsTransferConfiguration = null,
      virtualPrivateCloudConfiguration = null,
      batchRequestTimeoutConfiguration = null,
      referenceFileToDiskImageMappingOpt = None,
      dockerImageToCacheDiskImageMappingOpt = None,
      checkpointingInterval = 1 second)

    PipelinesApiBackendLifecycleActorFactory.robustBuildAttributes(() => attributes) shouldBe attributes
  }

  it should "retry retryable failures only" in {
    val retryMessage = "We encountered an internal error. Please try again."
    val fails = Table(
      ("function", "failing attempt", "description"),
      (() => throw new RuntimeException(), 1, "no exception message"),
      (() => throw new RuntimeException("non retryable failure"), 1, "not a retryable message"),
      (() => throw new RuntimeException(retryMessage), 3, "retryable message"),
      (() => throw new Error(retryMessage), 1, "error not exception")
    )
    forAll(fails) { (fn, failingAttempt, description) =>
      val e = the [RuntimeException] thrownBy {
        PipelinesApiBackendLifecycleActorFactory.robustBuildAttributes(fn)
      }
      withClue(description) {
        e.getMessage should startWith(s"Failed to build PipelinesApiConfigurationAttributes on attempt $failingAttempt of 3")
      }
    }
  }
}
