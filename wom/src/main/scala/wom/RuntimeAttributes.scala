package wom

import wom.expression.WomExpression

object RuntimeAttributesKeys {
  val DockerKey = "docker"
  val ContainerKey = "container" // New for WDL 1.1, preferred over "docker"
  val MaxRetriesKey = "maxRetries"

  val CpuKey = "cpu"
  val CpuPlatformKey = "cpuPlatform"

  val GpuKey = "gpuCount"
  val GpuTypeKey = "gpuType"

  val fuseMountKey = "fuseMount"
  val jobTimeoutKey = "jobTimeout"

  val MemoryKey = "memory"
  val FailOnStderrKey = "failOnStderr"
  val ContinueOnReturnCodeKey = "continueOnReturnCode"

  // WDL 1.1 required hints
  // Semantically, this is the same as continueOnReturnCode as the two attributes are combined at the parsing stage
  val ReturnCodesKey = "returnCodes"
  val GpuRequiredKey = "gpu"

  // WDL 1.1 reserved hints
  // "Reserved" means optional for an engine to implement, but it must follow the spec if it does
  // Cromwell supports `runtime.inputs.<qualified path>.localizationOptional` for single files
  // as well as `runtime.localizationOptional` to apply to all files in a task
  val LocalizationOptional = "localizationOptional"
  val Inputs = "inputs"
  val Outputs = "outputs"

  val sharedMemoryKey = "sharedMemorySize"
}

case class RuntimeAttributes(attributes: Map[String, WomExpression])
