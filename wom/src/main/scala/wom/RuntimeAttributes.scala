package wom

import wom.expression.WomExpression

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
}

case class RuntimeAttributes(attributes: Map[String, WomExpression])
