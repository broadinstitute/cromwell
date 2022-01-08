package cromwell.backend.impl.bcs

import org.scalatest.TryValues._

class BcsMountSpec extends BcsTestUtilSpec {
  behavior of s"BcsMountSpec"
  val ossObject = "oss://bcs-test/bcs-dir/"
  val localFile = "/home/admin/local-dir/"

  it should "be an input mount if src starts with oss://" in {
    var writeSupport = true
    var entryString = s"$ossObject $localFile $writeSupport"
    var entry = BcsMount.parse(entryString).success.value


    entry shouldBe a [BcsInputMount]
    BcsMount.toString(entry.src) shouldEqual ossObject
    BcsMount.toString(entry.dest) shouldEqual localFile
    entry.writeSupport shouldEqual writeSupport

    writeSupport = false

    entryString = s"$ossObject $localFile $writeSupport"
    entry = BcsMount.parse(entryString).success.value
    entry shouldBe a [BcsInputMount]
    BcsMount.toString(entry.src) shouldEqual ossObject
    BcsMount.toString(entry.dest) shouldEqual localFile
    entry.writeSupport shouldEqual writeSupport
  }

  it should "be an output mount if dest starts with oss://" in {
    var writeSupport = true
    var entryString = s"$localFile $ossObject $writeSupport"
    var entry = BcsMount.parse(entryString).success.value


    entry shouldBe a [BcsOutputMount]
    BcsMount.toString(entry.src) shouldEqual localFile
    BcsMount.toString(entry.dest) shouldEqual ossObject
    entry.writeSupport shouldEqual writeSupport

    writeSupport = false

    entryString = s"$localFile $ossObject $writeSupport"
    entry = BcsMount.parse(entryString).success.value
    entry shouldBe a [BcsOutputMount]
    BcsMount.toString(entry.src) shouldEqual localFile
    BcsMount.toString(entry.dest) shouldEqual ossObject
    entry.writeSupport shouldEqual writeSupport
  }

  it should "throw if src and dest are both local files or oss files" in {
    val writeSupport = true
    var entryString = s"$localFile $localFile $writeSupport"
    a [UnsupportedOperationException] should be thrownBy BcsMount.parse(entryString).failure.get

    entryString = s"$ossObject $ossObject $writeSupport"
    a [UnsupportedOperationException] should be thrownBy BcsMount.parse(entryString).failure.get
  }
}
