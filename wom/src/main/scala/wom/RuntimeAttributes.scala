package wom

import wom.expression.WomExpression

object RuntimeAttributesKeys {
  val DockerKey = "docker"
  /**
    * Equivalent to CPUMinKey
    */
  val CpuKey = "cpu"
  val CpuMinKey = "cpuMin"
  val CpuMaxKey = "cpuMax"
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
