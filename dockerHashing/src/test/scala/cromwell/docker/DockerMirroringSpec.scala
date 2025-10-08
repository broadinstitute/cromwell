package cromwell.docker

import com.typesafe.config.ConfigFactory
import cromwell.core.TestKitSuite
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks._

class DockerMirroringSpec extends TestKitSuite with AnyFlatSpecLike with Matchers {

  behavior of "DockerMirroring config parsing"

  it should "parse empty config correctly" in {
    val config = ConfigFactory.parseString("")
    DockerMirroring.fromConfig(config) shouldBe None
  }

  it should "parse dockerhub config correctly" in {
    val config = ConfigFactory.parseString("""
                                             |docker-mirror {
                                             |  dockerhub {
                                             |    enabled: true
                                             |    address: "foo.bar"
                                             |  }
                                             |}
                                             |""".stripMargin)
    DockerMirroring.fromConfig(config) shouldBe Some(DockerMirroring(mirrors = List(DockerHubMirror("foo.bar"))))
  }

  it should "parse no-address dockerhub config correctly" in {
    val config = ConfigFactory.parseString("""
                                             |docker-mirror {
                                             |  dockerhub {
                                             |    enabled: true
                                             |  }
                                             |}
                                             |""".stripMargin)
    DockerMirroring.fromConfig(config) shouldBe None
  }

  it should "parse disabled config correctly" in {
    val config = ConfigFactory.parseString("""
                                             |docker-mirror {
                                             |  dockerhub {
                                             |    enabled: false
                                             |    address: "foo.bar"
                                             |  }
                                             |}
                                             |""".stripMargin)
    DockerMirroring.fromConfig(config) shouldBe None
  }

  it should "parse non-dockerhub config (unsupported) correctly" in {
    val config = ConfigFactory.parseString("""
                                             |docker-mirror {
                                             |  quay {
                                             |    enabled: true
                                             |    address: "foo.bar"
                                             |  }
                                             |}
                                             |""".stripMargin)
    DockerMirroring.fromConfig(config) shouldBe None
  }

  behavior of "DockerMirroring"

  it should "apply a mirror with multiple path components" in {
    val mirroring = DockerMirroring(mirrors = List(DockerHubMirror("my.mirror.io/more/stuff")))
    val imgOpt = DockerImageIdentifier.fromString("docker.io/broadinstitute/cromwell")
    val img = imgOpt.get
    val m = mirroring.mirrorImage(img)
    m.map(_.fullName) shouldBe Some("my.mirror.io/more/stuff/broadinstitute/cromwell:latest")
    m.flatMap(_.host) shouldBe Some("my.mirror.io")
    m.flatMap(_.repository) shouldBe Some("more/stuff/broadinstitute")
    m.map(_.image) shouldBe Some("cromwell")
  }

  val mirroring = DockerMirroring(mirrors = List(DockerHubMirror("my.mirror.io")))

  val dockerMirrorInputs = Table(
    ("testName", "origImg", "mirrorResult"),
    ("mirror a Dockerhub image with default repository", "ubuntu:latest", Some("my.mirror.io/library/ubuntu:latest")),
    ("mirror a Dockerhub image with explicit repository",
     "broadinstitute/cromwell:v100",
     Some("my.mirror.io/broadinstitute/cromwell:v100")
    ),
    ("mirror a docker.io Dockerhub image",
     "docker.io/broadinstitute/cromwell",
     Some("my.mirror.io/broadinstitute/cromwell:latest")
    ),
    ("not mirror a GCR image", "gcr.io/broad-dsde-cromwell-dev/cromwell-drs-localizer", None)
  )

  forAll(dockerMirrorInputs) { (testName, origImg, mirrorResult) =>
    it should testName in {
      val imgOpt = DockerImageIdentifier.fromString(origImg)
      val img = imgOpt.get
      val m = mirroring.mirrorImage(img)
      val name = m.map(_.fullName)
      name shouldBe mirrorResult
    }
  }
}
