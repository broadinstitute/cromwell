package cromwell.core.callcaching.docker

import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}

class DockerImageIDSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {
  behavior of "DockerImageID"
  
  it should "parse valid docker images" in {
    val valid = Table(
      ("sourceString",                        "host",            "repo",          "image",     "reference"),
      
      // Without tags -> latest
      ("ubuntu",                      "registry-1.docker.io",   "library",        "ubuntu",     "latest"),
      ("broad/cromwell",              "registry-1.docker.io",   "broad",          "cromwell",   "latest"),
      ("index.docker.io/ubuntu",      "index.docker.io",        "library",        "ubuntu",     "latest"),
      ("broad/cromwell/submarine",    "registry-1.docker.io",   "broad/cromwell", "submarine",  "latest"),
      ("gcr.io/google/alpine",        "gcr.io",                 "google",         "alpine",     "latest"),

      // With tags
      ("ubuntu:latest",               "registry-1.docker.io",   "library",        "ubuntu",     "latest"),
      ("ubuntu:1235-SNAP",            "registry-1.docker.io",   "library",        "ubuntu",     "1235-SNAP"),
      ("ubuntu:V3.8-5_1",             "registry-1.docker.io",   "library",        "ubuntu",     "V3.8-5_1")
    )
    
    forAll(valid) { (dockerString, host, repo, image, reference) =>
      val imageId = DockerImageIdentifier.fromString(dockerString)
      imageId.isSuccess shouldBe true
      val successfulId = imageId.get
      successfulId.host shouldBe host
      successfulId.repository shouldBe repo
      successfulId.image shouldBe image
      successfulId.reference shouldBe reference
    }

    // With digest
    val withDigest = DockerImageIdentifier.fromString("ubuntu@sha256:45168651")
    withDigest.isSuccess shouldBe true
    val successfulDigest = withDigest.get
    successfulDigest.host shouldBe "registry-1.docker.io"
    successfulDigest.repository shouldBe "library"
    successfulDigest.image shouldBe "ubuntu"
    successfulDigest.reference shouldBe "sha256:45168651"
    successfulDigest.isInstanceOf[DockerImageIdentifierWithHash] shouldBe true
    successfulDigest.asInstanceOf[DockerImageIdentifierWithHash].hash shouldBe DockerHashResult("sha256", "45168651")
  }

  it should "not parse invalid docker images" in {
    val invalid = List(
      "_notvalid:latest",
      "NotValid:latest",
      "not:_valid",
      "/not:valid",
      "not%valid"
    )
    
    invalid foreach { image =>
      DockerImageIdentifier.fromString(image).isSuccess shouldBe false
    }
  }
}
