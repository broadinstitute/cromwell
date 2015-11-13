package cromwell.util.docker

import org.scalatest.FlatSpec

class DockerTestSpec extends FlatSpec {
  behavior of "DockerTest"

  it should "run" in {
    DockerTest.run()
  }
}
