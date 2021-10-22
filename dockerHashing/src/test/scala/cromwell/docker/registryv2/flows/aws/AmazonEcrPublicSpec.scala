package cromwell.docker.registryv2.flows.aws

import cromwell.core.TestKitSuite
import cromwell.docker.{DockerImageIdentifier, DockerRegistryConfig}
import org.scalatest.BeforeAndAfter
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.scalatestplus.mockito.MockitoSugar

class AmazonEcrPublicSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with MockitoSugar with BeforeAndAfter{
  val goodUri = "public.ecr.aws/amazonlinux/amazonlinux:latest"
  val badUri = "ubuntu:latest"
  val registry = new AmazonEcrPublic(DockerRegistryConfig.default)


  it should "Accept good URI" in {
    val dockerImageIdentifier = DockerImageIdentifier.fromString(goodUri).get
    registry.accepts(dockerImageIdentifier) shouldEqual(true)
  }

  it should "NOT accept bad URI" in {
    val dockerImageIdentifier = DockerImageIdentifier.fromString(badUri).get
    registry.accepts(dockerImageIdentifier) shouldEqual(false)
  }
}
