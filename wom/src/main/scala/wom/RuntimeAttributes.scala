package wom

import wom.expression.WomExpression
import wom.types.{WomMapType, WomStringType}

object RuntimeAttributesKeys {
  val DockerKey = "docker"
  val MaxRetriesKey = "maxRetries"

  val CpuKey = "cpu"
  val CpuPlatformKey = "cpuPlatform"

  val GpuKey = "gpuCount"
  val GpuTypeKey = "gpuType"

  val MemoryKey = "memory"
  val FailOnStderrKey = "failOnStderr"
  val ContinueOnReturnCodeKey = "continueOnReturnCode"

  object CloudProviders {
    val GcpKey = "gcp"
    val AzureKey = "azure"
  }
}

sealed trait CloudProvider {
  def runtimeKey: String
}

object CloudProvider {
  def all: Seq[CloudProvider] = Seq(GcpProvider, AzureProvider)
}

object GcpProvider extends CloudProvider {
  override def runtimeKey: String = "gcp"
}
object AzureProvider extends CloudProvider {
  override def runtimeKey: String = "azure"
}

case class RuntimeAttributes(attributes: Map[String, WomExpression]) {
  def forProvider(provider: CloudProvider): RuntimeAttributes = {
    attributes(provider.runtimeKey) match {
      case maybeAttributes: WomMap =>
        if (maybeAttributes.keyType == WomStringType)
          CloudProvider.all.
        else
          this
      case _ => {
        // If a user gives us "azure": "banana", "azure": 42, or anything else other than a map, ignore it
        this
      }
    }
  }
}
