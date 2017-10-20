package cromwell.core

import org.scalatest.Tag

object Tags {
  object DockerTest extends Tag("DockerTest")
  object IntegrationTest extends Tag("CromwellIntegrationTest")
  object DbmsTest extends Tag("DbmsTest")
  object PostWomTest extends Tag("PostWomTest")
}
