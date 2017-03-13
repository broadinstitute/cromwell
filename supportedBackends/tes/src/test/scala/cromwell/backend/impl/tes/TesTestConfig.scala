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
      |  # The keys below have been commented out as they are optional runtime attributes.
      |  # dockerWorkingDir
      |  # docker
      |}
      |""".stripMargin

  val backendConfig = ConfigFactory.parseString(backendConfigString)
}