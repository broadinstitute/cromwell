package cromwell.backend.google.batch.models

import cromwell.backend.google.batch.models.GpuResource.GpuType
import cromwell.backend.google.batch.models.GcpBatchTestConfig.gcpBatchConfiguration
import org.scalatest.matchers.should.Matchers
import org.scalatest.wordspec.AnyWordSpecLike
import wom.values.{WomFloat, WomInteger, WomSingleFile, WomString, WomValue}

class GcpBatchGpuAttributesSpec extends AnyWordSpecLike with Matchers with GcpBatchRuntimeAttributesSpecsMixin {

  val validGpuTypes = List(
    (Option(WomString("nvidia-tesla-k80")), Option(GpuType.NVIDIATeslaK80)),
    (Option(WomString("nvidia-tesla-p100")), Option(GpuType.NVIDIATeslaP100)),
    (Option(WomString("custom-gpu-24601")), Option(GpuType("custom-gpu-24601"))),
    (None, None)
  )
  val invalidGpuTypes = List(WomSingleFile("nvidia-tesla-k80"), WomInteger(100))

  val validGpuCounts = List(
    (Option(WomInteger(1)), Option(1)),
    (Option(WomInteger(100)), Option(100)),
    (None, None)
  )
  val invalidGpuCounts = List(WomString("ten"), WomFloat(1.0))

  validGpuTypes foreach { case (validGpuType, expectedGpuTypeValue) =>
    validGpuCounts foreach { case (validGpuCount, expectedGpuCountValue) =>
      s"validate the valid gpu type '$validGpuType' and count '$validGpuCount'" in {
        val runtimeAttributes = Map(
          "docker" -> WomString("ubuntu:latest")
        ) ++ validGpuType.map(t => "gpuType" -> t) ++ validGpuCount.map(c => "gpuCount" -> c)

        val actualRuntimeAttributes =
          toBatchRuntimeAttributes(runtimeAttributes, emptyWorkflowOptions, gcpBatchConfiguration)

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

        assertBatchRuntimeAttributesFailedCreation(runtimeAttributes,
                                                   s"Invalid gpu count. Expected positive Int but got"
        )
      }
    }
  }

  invalidGpuTypes foreach { invalidGpuType =>
    invalidGpuCounts foreach { invalidGpuCount =>
      s"not validate a invalid gpu type '$invalidGpuType' and invalid gpu count '$invalidGpuCount'" in {
        val runtimeAttributes: Map[String, WomValue] = Map(
          "docker" -> WomString("ubuntu:latest")
        ) + ("gpuType" -> invalidGpuType) + ("gpuCount" -> invalidGpuCount)

        assertBatchRuntimeAttributesFailedCreation(runtimeAttributes,
                                                   s"Invalid gpu count. Expected positive Int but got"
        )
      }
    }
  }
}
