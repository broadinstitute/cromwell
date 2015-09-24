package cromwell

import org.scalatest.Tag

object CromwellSpec {
  val BackendType = cromwell.parser.BackendType.LOCAL

  object DockerTest extends Tag("DockerTest")
}
