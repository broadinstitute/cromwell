package cromwell

import com.typesafe.config.ConfigFactory

object CromwellSpec {
  val BackendConfText =
    """
      |backend {
      |   backend = "local"
      |}
    """.stripMargin
  val Config = ConfigFactory.parseString(CromwellSpec.BackendConfText).getConfig("backend")
}
