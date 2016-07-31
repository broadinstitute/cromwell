package cromwell.backend.impl.sfs.config

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
  val CheckAliveConfig = "check-alive"
  val RuntimeAttributesConfig = "runtime-attributes"
  val JobIdRegexConfig = "job-id-regex"
  val RunInBackgroundConfig = "run-in-background"

  /*
  Runtime attributes that may be specified within the RuntimeAttributesConfig.
   */
  val DockerRuntimeAttribute = "docker"
  val CpuRuntimeAttribute = "cpu"
  val MemoryRuntimeAttribute = "memory"
  // See: MemoryDeclarationValidation
  val MemoryRuntimeAttributePrefix = "memory_"

  /*
  List of task names used internally.
  NOTE: underscore separated
   */
  val RuntimeAttributesTask = "runtime_attributes"
  val SubmitTask = "submit"
  val SubmitDockerTask = "submit_docker"
  val KillTask = "kill"
  val CheckAliveTask = "check_alive"

  /*
  Inputs passed into the submit and kill commands.
  NOTE: underscore separated.
   */
  val JobNameInput = "job_name"
  val CwdInput = "cwd"
  val DockerCwdInput = "docker_cwd"
  val StdoutInput = "out"
  val StderrInput = "err"
  val ScriptInput = "script"
  val JobIdInput = "job_id"
}
