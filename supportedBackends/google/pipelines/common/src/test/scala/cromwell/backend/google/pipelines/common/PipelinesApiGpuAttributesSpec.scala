package cromwell.backend.google.pipelines.common

import cromwell.backend.google.pipelines.common.GpuResource.GpuType
import cromwell.backend.google.pipelines.common.PipelinesApiTestConfig.papiConfiguration
import eu.timepit.refined.refineMV
import org.scalatest.{Matchers, WordSpecLike}
import wom.values.{WomFloat, WomInteger, WomSingleFile, WomString, WomValue}

class PipelinesApiGpuAttributesSpec
  extends WordSpecLike
    with Matchers
    with PipelinesApiRuntimeAttributesSpecsMixin {

  val validGpuTypes = List(
    (Option(WomString("nvidia-tesla-k80")), Option(GpuType.NVIDIATeslaK80)),
    (Option(WomString("nvidia-tesla-p100")), Option( GpuType.NVIDIATeslaP100)),
    (Option(WomString("custom-gpu-24601")), Option( GpuType("custom-gpu-24601"))),
    (None, None))
  val invalidGpuTypes = List(
    WomSingleFile("nvidia-tesla-k80"),
    WomInteger(100))

  val validGpuCounts = List(
    (Option(WomInteger(1)), Option(1)),
    (Option(WomInteger(100)), Option(100)),
    (None, None)
  )
  val invalidGpuCounts = List(
    WomString("ten"),
    WomFloat(1.0))

  validGpuTypes foreach { case (validGpuType, expectedGpuTypeValue) =>
    validGpuCounts foreach { case (validGpuCount, expectedGpuCountValue) =>
      s"validate the valid gpu type '$validGpuType' and count '$validGpuCount'" in {
        val runtimeAttributes = Map(
          "docker" -> WomString("ubuntu:latest")
        ) ++ validGpuType.map(t => "gpuType" -> t) ++ validGpuCount.map(c => "gpuCount" -> c)

        val actualRuntimeAttributes = toPapiRuntimeAttributes(runtimeAttributes, emptyWorkflowOptions, papiConfiguration)

        expectedGpuTypeValue match {
          case Some(v) => actualRuntimeAttributes.gpuResource.exists(_.gpuType == v)
          case None => actualRuntimeAttributes.gpuResource.foreach(_.gpuType == GpuType.DefaultGpuType)
        }

        expectedGpuCountValue match {
          case Some(v) => actualRuntimeAttributes.gpuResource.exists(_.gpuCount.value == v)
          case None => actualRuntimeAttributes.gpuResource.foreach(_.gpuCount.value == GpuType.DefaultGpuCount.value)
        }

      }
    }

    invalidGpuCounts foreach { invalidGpuCount =>
      s"not validate a valid gpu type '$validGpuType' but an invalid gpu count '$invalidGpuCount'" in {
        val runtimeAttributes: Map[String, WomValue] = Map(
          "docker" -> WomString("ubuntu:latest")
        ) ++ validGpuType.map(t => "gpuType" -> t) + ("gpuCount" -> invalidGpuCount)

        assertJesRuntimeAttributesFailedCreation(
          runtimeAttributes,
          s"Invalid gpu count. Expected positive Int but got")
      }
    }
  }

  invalidGpuTypes foreach { invalidGpuType =>
    invalidGpuCounts foreach { invalidGpuCount =>
      s"not validate a invalid gpu type '$invalidGpuType' and invalid gpu count '$invalidGpuCount'" in {
        val runtimeAttributes: Map[String, WomValue] = Map(
          "docker" -> WomString("ubuntu:latest")
        ) + ("gpuType" -> invalidGpuType) + ("gpuCount" -> invalidGpuCount)

        assertJesRuntimeAttributesFailedCreation(
          runtimeAttributes,
          s"Invalid gpu count. Expected positive Int but got")
      }
    }
  }

  "validate a valid gpu entry (1)" in {
    val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "gpuCount" -> WomInteger(1), "gpuType" -> WomString("nvidia-tesla-k80"))
    val expectedRuntimeAttributes = expectedDefaults.copy(gpuResource = Option(GpuResource(gpuCount = refineMV(1), gpuType = GpuType.NVIDIATeslaK80)))
    assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
  }

  "validate a valid gpu entry (2)" in {
    val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "gpuCount" -> WomInteger(2), "gpuType" -> WomString("nvidia-tesla-p100"))
    val expectedRuntimeAttributes = expectedDefaults.copy(gpuResource = Option(GpuResource(gpuCount = refineMV(2), gpuType = GpuType.NVIDIATeslaP100)))
    assertJesRuntimeAttributesSuccessfulCreation(runtimeAttributes, expectedRuntimeAttributes)
  }

  // Missing gpu type
  "fail to validate an invalid gpu entry (1)" in {
    val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "gpuCount" -> WomInteger(1))
    assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Please specify a GPU type: nvidia-tesla-p100, nvidia-tesla-k80")
  }

  // Missing gpu count
  "fail to validate an invalid gpu entry (2)" in {
    val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "gpuType" -> WomString("nvidia-tesla-p100"))
    assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Please specify how many GPU should be attached to the instance.")
  }

  // unrecoginzed gpu type
  "fail to validate an invalid gpu entry (3)" in {
    val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "gpuCount" -> WomInteger(1), "gpuType" -> WomString("not-a-gpu"))
    assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "not-a-gpu is not a supported GPU type. Supported types are nvidia-tesla-k80, nvidia-tesla-p100")
  }

  // gpu count is not an int
  "fail to validate an invalid gpu entry (4)" in {
    val runtimeAttributes = Map("docker" -> WomString("ubuntu:latest"), "gpuCount" -> WomString("value"))
    assertJesRuntimeAttributesFailedCreation(runtimeAttributes, "Expecting gpuCount runtime attribute to be an Integer")
  }


}
