package cromwell.backend.validation

import wom.values.WomValue

object RuntimeAttributesKeys {
  val DockerKey = "docker"
  val FailOnStderrKey = "failOnStderr"
  val ContinueOnReturnCodeKey = "continueOnReturnCode"
  val CpuKey = "cpu"
  val MemoryKey = "memory"
}

case class RuntimeKey[T](key: String) {
  def apply(attributes: Map[String, WomValue]):
}
