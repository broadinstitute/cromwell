package cromwell.backend.impl.local

import com.typesafe.config.ConfigFactory

trait HasLocalBackendConfig {
  protected lazy val localBackendConfig = ConfigFactory.parseString(
    """
      |{
      |  root: "local-cromwell-executions"
      |  filesystems {
      |    local {
      |      localization: [
      |        "hard-link", "soft-link", "copy"
      |      ]
      |    }
      |  }
      |}
    """.stripMargin)
}
