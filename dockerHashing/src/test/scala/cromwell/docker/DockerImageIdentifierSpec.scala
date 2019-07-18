package cromwell.docker

import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}

class DockerImageIdentifierSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {
  behavior of "DockerImageID"

  it should "parse valid docker images" in {
    val valid = Table(
      ("sourceString",                            "host",             "repo",                 "image",     "reference"),
      // Without tags -> latest
      ("ubuntu",                                  None,               None,                   "ubuntu",     "latest"),
      ("broad/cromwell",                          None,               Option("broad"),          "cromwell",   "latest"),
      ("index.docker.io/ubuntu",         Option("index.docker.io"),   None,                   "ubuntu",     "latest"),
      ("broad/cromwell/submarine",                None,               Option("broad/cromwell"), "submarine",  "latest"),
      ("gcr.io/google/slim",              Option("gcr.io"),         Option("google"),         "slim",     "latest"),
      // With tags
      ("ubuntu:latest",                           None,               None,                   "ubuntu",     "latest"),
      ("ubuntu:1235-SNAP",                        None,               None,                   "ubuntu",     "1235-SNAP"),
      ("ubuntu:V3.8-5_1",                         None,               None,                   "ubuntu",     "V3.8-5_1"),
      ("index.docker.io:9999/ubuntu:170904",  Option("index.docker.io:9999"), None,             "ubuntu",    "170904"),
      ("localhost:5000/capture/transwf:170904", Option("localhost:5000"), Option("capture"),      "transwf",    "170904"),
      ("quay.io/biocontainers/platypus-variant:0.8.1.1--htslib1.5_0", Option("quay.io"), Option("biocontainers"), "platypus-variant", "0.8.1.1--htslib1.5_0"),
      ("registry.cn-shanghai.aliyuncs.com/batchcompute/ubuntu:0.2", Option("registry.cn-shanghai.aliyuncs.com"), Option("batchcompute"), "ubuntu", "0.2")
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
    successfulDigest.repository shouldBe None
    successfulDigest.image shouldBe "ubuntu"
    successfulDigest.reference shouldBe "sha256:45168651"
    successfulDigest.isInstanceOf[DockerImageIdentifierWithHash] shouldBe true
    successfulDigest.asInstanceOf[DockerImageIdentifierWithHash].hash shouldBe DockerHashResult("sha256", "45168651")

    // With tag + digest
    val withTagDigest = DockerImageIdentifier.fromString("ubuntu:latest@sha256:45168651")
    withTagDigest.isSuccess shouldBe true
    val successfulTagDigest = withTagDigest.get
    successfulTagDigest.host shouldBe None
    successfulTagDigest.repository shouldBe None
    successfulTagDigest.image shouldBe "ubuntu"
    successfulTagDigest.reference shouldBe "latest@sha256:45168651"
    successfulTagDigest.isInstanceOf[DockerImageIdentifierWithHash] shouldBe true
    successfulTagDigest.asInstanceOf[DockerImageIdentifierWithHash].hash shouldBe DockerHashResult("sha256", "45168651")
  }

  it should "not parse invalid docker images" in {
    val invalid = List(
      "_notvalid:latest",
      "NotValid:latest",
      "not:_valid",
      "/not:valid",
      "not%valid",
      "not@sha256:digest:tag",
      "not:sha256:digest@tag"
    )

    invalid foreach { image =>
      DockerImageIdentifier.fromString(image).isSuccess shouldBe false
    }
  }
}
