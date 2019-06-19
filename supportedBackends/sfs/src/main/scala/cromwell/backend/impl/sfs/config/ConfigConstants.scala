package cromwell.backend.impl.sfs.config

import wom.RuntimeAttributesKeys

/**
  * A consolidated list of constants, a one-stop-shop for backend configuration to know what constants are floating
  * around.
  */
object ConfigConstants {
  /*
  List of config keys that may be used in the application.conf.
  NOTE: hyphen separated
   */
  val SubmitConfig = "submit"
  val SubmitDockerConfig = "submit-docker"
  val KillConfig = "kill"
  val KillDockerConfig = "kill-docker"
  val CheckAliveConfig = "check-alive"
  val ExitCodeTimeoutConfig = "exit-code-timeout-seconds"
  val RuntimeAttributesConfig = "runtime-attributes"
  val JobIdRegexConfig = "job-id-regex"
  val RunInBackgroundConfig = "run-in-background"

  /*
  Runtime attributes that may be specified within the RuntimeAttributesConfig.
   */
  val DockerRuntimeAttribute = RuntimeAttributesKeys.DockerKey
  val CpuRuntimeAttribute = RuntimeAttributesKeys.CpuKey
  val MemoryRuntimeAttribute = RuntimeAttributesKeys.MemoryKey
  val MemoryMinRuntimeAttribute = RuntimeAttributesKeys.MemoryMinKey
  val MemoryMaxRuntimeAttribute = RuntimeAttributesKeys.MemoryMaxKey
  // See: MemoryDeclarationValidation
  val MemoryRuntimeAttributePrefix = "memory_"
  val MemoryMinRuntimeAttributePrefix = "memoryMin_"
  val MemoryMaxRuntimeAttributePrefix = "memoryMax_"
  val DiskRuntimeAttribute = "disk"
  val DiskRuntimeAttributePrefix = "disk_"
  /*
  List of task names used internally.
  NOTE: underscore separated
   */
  val RuntimeAttributesTask = "runtime_attributes"
  val SubmitTask = "submit"
  val SubmitDockerTask = "submit_docker"
  val KillTask = "kill"
  val KillDockerTask = "kill_docker"
  val CheckAliveTask = "check_alive"

  /*
  Inputs passed into the submit and kill commands.
  NOTE: underscore separated.
   */
  val JobNameInput = "job_name"
  val CwdInput = "cwd"
  val DockerCwdInput = "docker_cwd"
  val DockerCidInput = "docker_cid"
  val DockerScriptInput = "docker_script"
  val DockerStdoutInput = "docker_out"
  val DockerStderrInput = "docker_err"
  val StdoutInput = "out"
  val StderrInput = "err"
  val ScriptInput = "script"
  val JobIdInput = "job_id"
  val JobShellInput = "job_shell"
}
