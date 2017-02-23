package cromwell.core.callcaching.docker

import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}

class DockerImageIdentifierSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {
  behavior of "DockerImageID"
  
  it should "parse valid docker images" in {
    val valid = Table(
      ("sourceString",                            "host",              "repo",          "image",     "reference"),
      
      // Without tags -> latest
      ("ubuntu",                                  None,               "library",        "ubuntu",     "latest"),
      ("broad/cromwell",                          None,               "broad",          "cromwell",   "latest"),
      ("index.docker.io/ubuntu",         Option("index.docker.io"),   "library",        "ubuntu",     "latest"),
      ("broad/cromwell/submarine",                None,               "broad/cromwell", "submarine",  "latest"),
      ("gcr.io/google/alpine",              Option("gcr.io"),         "google",         "alpine",     "latest"),

      // With tags
      ("ubuntu:latest",                           None,               "library",        "ubuntu",     "latest"),
      ("ubuntu:1235-SNAP",                        None,               "library",        "ubuntu",     "1235-SNAP"),
      ("ubuntu:V3.8-5_1",                         None,               "library",        "ubuntu",     "V3.8-5_1")
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
    successfulDigest.host shouldBe None
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
