package wom

import wom.expression.WomExpression

object RuntimeAttributesKeys {
  val DockerKey = "docker"
  val MaxRetriesKey = "maxRetries"
  /**
    * Equivalent to CPUMinKey
    */
  val CpuKey = "cpu"
  val CpuPlatformKey = "cpuPlatform"
  val CpuMinKey = "cpuMin"
  val CpuMaxKey = "cpuMax"
  /**
    * Equivalent to GPUMinKey
    */
  val GpuKey = "gpuCount"
  val GpuMinKey = "gpuCountMin"
  val GpuMaxKey = "gpuCountMax"
  val GpuTypeKey = "gpuType"
  val DnaNexusInputDirMinKey = "dnaNexusInputDirMin"
  /**
    * Equivalent to MemoryMinKey
    */
  val MemoryKey = "memory"
  val MemoryMinKey = "memoryMin"
  val MemoryMaxKey = "memoryMax"
  val TmpDirMinKey = "tmpDirMin"
  val TmpDirMaxKey = "tmpDirMax"
  val OutDirMinKey = "outDirMin"
  val OutDirMaxKey = "outDirMax"
  val FailOnStderrKey = "failOnStderr"
  val ContinueOnReturnCodeKey = "continueOnReturnCode"
}

case class RuntimeAttributes(attributes: Map[String, WomExpression])
