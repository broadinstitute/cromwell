package cromwell.backend.google.pipelines.common

import eu.timepit.refined.numeric.Positive
import eu.timepit.refined.refineV
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

import scala.concurrent.duration._
import scala.language.postfixOps

class PipelinesApiBackendLifecycleActorFactorySpec
    extends AnyFlatSpecLike
    with Matchers
    with TableDrivenPropertyChecks {

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
      qps = refineV[Positive](1).toOption.get,
      cacheHitDuplicationStrategy = null,
      requestWorkers = refineV[Positive](1).toOption.get,
      pipelineTimeout = 1 second,
      logFlushPeriod = Option(1 second),
      gcsTransferConfiguration = null,
      virtualPrivateCloudConfiguration = null,
      batchRequestTimeoutConfiguration = null,
      referenceFileToDiskImageMappingOpt = None,
      dockerImageToCacheDiskImageMappingOpt = None,
      checkpointingInterval = 1 second
    )

    PipelinesApiBackendLifecycleActorFactory.robustBuildAttributes(() => attributes) shouldBe attributes
  }

  {
    // The message string actually observed during construction failure.
    val actualRetryMessage = s"We encountered an internal error. Please try again."
    val fails = Table(
      ("attempts", "description", "function"),
      (3, "no exception message", () => throw new RuntimeException()),
      (3, "exception message", () => throw new RuntimeException(actualRetryMessage)),
      (1, "error not exception, no message", () => throw new Error()),
      (1, "error not exception", () => throw new Error(actualRetryMessage))
    )
    forAll(fails) { (attempts, description, function) =>
      it should s"$description: make $attempts attribute creation attempts before giving up" in {
        val e = the[RuntimeException] thrownBy {
          PipelinesApiBackendLifecycleActorFactory.robustBuildAttributes(function,
                                                                         initialIntervalMillis = 1,
                                                                         maxIntervalMillis = 5
          )
        }
        e.getMessage should startWith(s"Failed to build PipelinesApiConfigurationAttributes on attempt $attempts of 3")
      }
    }
  }
}
