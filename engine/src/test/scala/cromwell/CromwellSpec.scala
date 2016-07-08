package cromwell

import com.typesafe.config.ConfigFactory
import org.scalatest.Tag

object CromwellSpec {
  val BackendConfText =
    """
      |backend {
      |   backend = "local"
      |}
    """.stripMargin
  val Config = ConfigFactory.parseString(CromwellSpec.BackendConfText).getConfig("backend")

  object DockerTest extends Tag("DockerTest")
  object PostMVP extends Tag("PostMVP")
  object IntegrationTest extends Tag("CromwellIntegrationTest")
  object DbmsTest extends Tag("DbmsTest")
  object MainSpecWdlRunner extends Tag("MainSpecWdlRunner")
}
