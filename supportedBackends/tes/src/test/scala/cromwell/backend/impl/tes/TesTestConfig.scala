package cromwell.backend.impl.tes

import com.typesafe.config.ConfigFactory

object TesTestConfig {

  private val backendConfigString =
    """
      |root = "local-cromwell-executions"
      |dockerRoot = "/cromwell-executions"
      |endpoint = "http://127.0.0.1:9000/v1/jobs"
      |
      |default-runtime-attributes {
      |  cpu: 1
      |  failOnStderr: false
      |  continueOnReturnCode: 0
      |  memory: "2 GB"
      |  disk: "2 GB"
      |  preemptible: false
      |  # The keys below have been commented out as they are optional runtime attributes.
      |  # dockerWorkingDir
      |  # docker
      |}
      |""".stripMargin

  val backendConfig = ConfigFactory.parseString(backendConfigString)

  private val backendConfigStringWithBackendParams =
    """
      |root = "local-cromwell-executions"
      |dockerRoot = "/cromwell-executions"
      |endpoint = "http://127.0.0.1:9000/v1/jobs"
      |use_tes_11_preview_backend_parameters = true
      |
      |default-runtime-attributes {
      |  cpu: 1
      |  failOnStderr: false
      |  continueOnReturnCode: 0
      |  memory: "2 GB"
      |  disk: "2 GB"
      |  preemptible: false
      |  # The keys below have been commented out as they are optional runtime attributes.
      |  # dockerWorkingDir
      |  # docker
      |}
      |""".stripMargin

  val backendConfigWithBackendParams = ConfigFactory.parseString(backendConfigStringWithBackendParams)

  private val backendConfigStringWithBackoffs =
    """
      |root = "local-cromwell-executions"
      |dockerRoot = "/cromwell-executions"
      |endpoint = "http://127.0.0.1:9000/v1/jobs"
      |
      |default-runtime-attributes {
      |  cpu: 1
      |  failOnStderr: false
      |  continueOnReturnCode: 0
      |  memory: "2 GB"
      |  disk: "2 GB"
      |  preemptible: false
      |  # The keys below have been commented out as they are optional runtime attributes.
      |  # dockerWorkingDir
      |  # docker
      |}
      |
      |poll-backoff {
      |  min: "5 seconds"
      |  max: "1 minute"
      |  multiplier: 2.5
      |  randomization-factor: .7
      |}
      |
      |execute-or-recover-backoff {
      |  min: "3 minutes"
      |  max: "1 hours"
      |  multiplier: 5
      |  randomization-factor: .1
      |}
      |
      |""".stripMargin

  val backendConfigWithBackoffs = ConfigFactory.parseString(backendConfigStringWithBackoffs)

  private val backendConfigStringWithInvalidBackoffs =
    """
      |root = "local-cromwell-executions"
      |dockerRoot = "/cromwell-executions"
      |endpoint = "http://127.0.0.1:9000/v1/jobs"
      |
      |default-runtime-attributes {
      |  cpu: 1
      |  failOnStderr: false
      |  continueOnReturnCode: 0
      |  memory: "2 GB"
      |  disk: "2 GB"
      |  preemptible: false
      |  # The keys below have been commented out as they are optional runtime attributes.
      |  # dockerWorkingDir
      |  # docker
      |}
      |
      |poll-backoff {
      |  min: "5 seconds"
      |  max: "1 minute"
      |  multiplier: 2.5
      |  # missing randomization-factor
      |}
      |
      |""".stripMargin

  val backendConfigWithInvalidBackoffs = ConfigFactory.parseString(backendConfigStringWithInvalidBackoffs)
}
