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

  // New for WDL 1.1
  // Semantically, this is the same as continueOnReturnCode as the two attributes are combined at the parsing stage
  val ReturnCodesKey = "returnCodes"
}

case class RuntimeAttributes(attributes: Map[String, WomExpression])
