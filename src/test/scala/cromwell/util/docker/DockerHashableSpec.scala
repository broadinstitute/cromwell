package cromwell.util.docker

import cromwell.util.AggregatedException
import org.scalatest.{Matchers, FlatSpec}
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table

class DockerHashableSpec extends FlatSpec with Matchers {
  behavior of "DockerHashable"

  it should "create hashes for valid hashables" in {
    val hashables = Table(
      ("hashable",
        "hashType",
        "hashString")
      , (DockerRegistryImageId("1111222233334444555566667777888899990000aaaabbbbccccddddeeeeffff"),
        "imageId-sha256",
        "1111222233334444555566667777888899990000aaaabbbbccccddddeeeeffff")
      , (DockerHubImageId(Seq(DockerHubPartialLayerId("12345678"))),
        "layerIds-sha256Part",
        "12345678")
      , (DockerHubImageId(Seq(DockerHubPartialLayerId("12345678"), DockerHubPartialLayerId("87654321"))),
        "layerIds-sha256Part-md5",
        "461595e4bdb090ce41e7818287954d86")
      , (DockerManifest(Seq(DockerFsLayer("sha256:1111222233334444555566667777888899990000aaaabbbbccccddddeeeeffff"))),
        "layerBlobs-sha256",
        "1111222233334444555566667777888899990000aaaabbbbccccddddeeeeffff")
      , (DockerManifest(Seq(
        DockerFsLayer("sha256:1111222233334444555566667777888899990000aaaabbbbccccddddeeeeffff"),
        DockerFsLayer("sha256:ffffeeeeddddccccbbbbaaaa0000999988887777666655554444333322221111"))),
        "layerBlobs-sha256-md5",
        "c5a37bb0fdc0838a15914686cb85c7b9")
      , (DockerDigestHashable("sha256:1111222233334444555566667777888899990000aaaabbbbccccddddeeeeffff"),
        "digest-sha256",
        "1111222233334444555566667777888899990000aaaabbbbccccddddeeeeffff")
    )

    forAll(hashables) { (hashable, hashType, hashString) =>
      val dockerHash = hashable.dockerHash.get
      dockerHash.hashType should be(hashType)
      dockerHash.hashString should be(hashString)
    }
  }

  it should "not create hashes for invalid hashables" in {
    val hashables = Table(
      ("hashable", "exceptionType", "exceptionMessage"),
      (DockerRegistryImageId("bad_hash"), an[AggregatedException[_]],
        "hashString 'bad_hash' is not valid: contains illegal character for hexBinary: bad_hash"),
      (DockerHubImageId(Seq.empty), an[IllegalArgumentException], "docker hashes is empty"),
      (DockerManifest(Seq.empty), an[IllegalArgumentException], "docker hashes is empty"),
      (DockerDigestHashable("sha256:bad_hash"), an[AggregatedException[_]],
        "hashString 'bad_hash' is not valid: contains illegal character for hexBinary: bad_hash")
    )

    forAll(hashables) { (hashable, exceptionType, exceptionMessage) =>
      val exception = hashable.dockerHash.failed.get
      exception should be(exceptionType)
      exception.getMessage should be(exceptionMessage)
    }
  }
}
