package wom.values

import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks
import wom.WomExpressionException
import wom.types._

import scala.util.{Success, Try}

class WomFileSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "WomFile"

  lazy val singleDir = WomUnlistedDirectory("single/dir")
  lazy val singleFile = WomSingleFile("single/file")
  lazy val globFile = WomGlobFile("glob/*")
  lazy val listedDir1 = WomUnlistedDirectory("listed/dir1")
  lazy val listedDir2 = WomUnlistedDirectory("listed/dir2")
  lazy val secondaryFile1 = WomSingleFile("secondary/file1")
  lazy val secondaryFile2 = WomSingleFile("secondary/file2")

  val mapFileTests = Table(
    ("description", "womFile", "expected"),
    ("a single directory", singleDir, WomUnlistedDirectory("prepend/single/dir")),
    ("a single file", singleFile, WomSingleFile("prepend/single/file")),
    ("a glob file", globFile, WomGlobFile("prepend/glob/*"))
  )

  forAll(mapFileTests) { (description, womFile, expected) =>
    it should s"map $description" in {
      womFile.mapFile("prepend/" + _) should be(expected)
    }
  }

  forAll(mapFileTests) { (description, womFile, expected) =>
    it should s"map $description with mapWomFile" in {
      womFile.mapWomFile("prepend/" + _.value) should be(expected)
    }
  }

  val flattenFileTests = Table(
    ("description", "womFile", "expected"),
    ("a single directory", singleDir, List(WomUnlistedDirectory("single/dir"))),
    ("a single file", singleFile, List(WomSingleFile("single/file"))),
    ("a glob file", globFile, List(WomGlobFile("glob/*")))
  )

  forAll(flattenFileTests) { (description, womFile, expected) =>
    it should s"flatten $description" in {
      womFile.flattenFiles should be(expected)
    }
  }

  it should "test valueString" in {
    singleDir.valueString should be("single/dir")
    singleFile.valueString should be("single/file")
    globFile.valueString should be("glob/*")
  }

  it should "test womString" in {
    singleDir.toWomString should be(""""single/dir"""")
    singleFile.toWomString should be(""""single/file"""")
    globFile.toWomString should be("""glob("glob/*")""")
  }

  val addTests = Table(
    ("description", "womFile", "expectedPrefix", "expectedSuffix"),
    ("a single directory", singleDir, WomString("prefix/single/dir"), WomUnlistedDirectory("single/dir/suffix")),
    ("a single file", singleFile, WomString("prefix/single/file"), WomSingleFile("single/file/suffix")),
    ("a glob file", globFile, WomString("prefix/glob/*"), WomGlobFile("glob/*/suffix"))
  )

  forAll(addTests) { (description, womFile, expectedPrefix, expectedSuffix) =>
    it should s"add a string prefix to $description" in {
      WomString("prefix/").add(womFile) should be(Success(expectedPrefix))
    }

    it should s"add a string suffix to $description" in {
      womFile.add(WomString("/suffix")) should be(Success(expectedSuffix))
    }

    it should s"add an optional string prefix to $description" in {
      WomOptionalValue(WomStringType, Option(WomString("prefix/"))).add(
        WomOptionalValue(womFile.womType, Option(womFile))
      ) should be(Success(expectedPrefix))
    }

    it should s"add an optional string suffix to $description" in {
      WomOptionalValue(womFile.womType, Option(womFile)).add(
        WomOptionalValue(WomStringType, Option(WomString("/suffix")))
      ) should be(Success(expectedSuffix))
    }

    it should s"add a string suffix to $description, trimming that suffix before addition" in {
      womFile.add(WomString("  \n\r\t /suffix  \n\r\t ")) should be(Success(expectedSuffix))
    }

    it should s"add an optional string suffix to $description, trimming that suffix before addition" in {
      WomOptionalValue(womFile.womType, Option(womFile)).add(
        WomOptionalValue(WomStringType, Option(WomString("  \n\r\t /suffix  \n\r\t ")))
      ) should be(Success(expectedSuffix))
    }

  }

   val addInvalidIntegerTests = Table(
    ("description", "womFile", "expected"),
    ("a single directory", singleDir, "Cannot perform operation: single/dir + WomInteger(42)"),
    ("a single file", singleFile, "Cannot perform operation: single/file + WomInteger(42)"),
    ("a glob file", globFile, "Cannot perform operation: glob/* + WomInteger(42)")
  )

  forAll(addInvalidIntegerTests) { (description, womFile, expected) =>
    it should s"fail to add an int to $description" in {
      womFile.add(WomInteger(42)).failed.get should have message expected
    }

    it should s"fail to add an optional int to $description" in {
      WomOptionalValue(womFile.womType, Option(womFile))
        .add(
          WomOptionalValue(WomIntegerType, Option(WomInteger(42)))
        )
        .failed
        .get should have message expected
    }
  }

  val womFileEqualsTests = Table(
    ("description", "womFileA", "womFileB", "expected"),
    ("a single directory matched to a similar directory", singleDir, WomUnlistedDirectory(singleDir.value), true),
    ("a single file matched to a similar file", singleFile, WomSingleFile(singleFile.value), true),
    ("a glob file matched to a similar glob", globFile, WomGlobFile(globFile.value), true),
    ("a single directory matched to a dissimilar directory",
     singleDir,
     WomUnlistedDirectory("should/not/match"),
     false
    ),
    ("a single file matched to a dissimilar file", singleFile, WomSingleFile("should/not/match"), false),
    ("a glob file matched to a dissimilar glob", globFile, WomGlobFile("should/not/match"), false)
  )

  forAll(womFileEqualsTests) { (description, womFileA, womFileB, expected) =>
    it should s"expect $expected comparing $description" in {
      womFileA.equals(womFileB) should be(Success(WomBoolean(expected)))
    }

    it should s"expect $expected symmetrically comparing $description" in {
      womFileB.equals(womFileA) should be(Success(WomBoolean(expected)))
    }

    it should s"expect $expected comparing (optionally) $description" in {
      WomOptionalValue(womFileA.womType, Option(womFileA)).equals(
        WomOptionalValue(womFileB.womType, Option(womFileB))
      ) should be(Success(WomBoolean(expected)))
    }

    it should s"expect $expected symmetrically comparing (optionally) $description" in {
      WomOptionalValue(womFileB.womType, Option(womFileB)).equals(
        WomOptionalValue(womFileA.womType, Option(womFileA))
      ) should be(Success(WomBoolean(expected)))
    }
  }

  val womFileEqualsFailures = Table(
    ("description", "womFileA", "womFileB"),
    ("a single file matched to a similar directory", singleFile, WomUnlistedDirectory(singleDir.value)),
    ("a glob file matched to a similar directory", globFile, WomUnlistedDirectory(singleDir.value)),
    ("a single directory matched to a similar file", singleDir, WomSingleFile(singleFile.value)),
    ("a glob file matched to a similar file", globFile, WomSingleFile(singleFile.value)),
    ("a single directory matched to a similar glob", singleDir, WomGlobFile(globFile.value)),
    ("a single file matched to a similar glob", singleFile, WomGlobFile(globFile.value))
  )

  forAll(womFileEqualsFailures) { (description, womFileA, womFileB) =>
    it should s"fail comparing $description" in {
      a[WomExpressionException] should be thrownBy {
        womFileA.equals(womFileB).get
      }
    }

    it should s"fail symmetrically comparing $description" in {
      a[WomExpressionException] should be thrownBy {
        womFileB.equals(womFileA).get
      }
    }

    it should s"fail comparing (optionally) $description" in {
      a[WomExpressionException] should be thrownBy {
        WomOptionalValue(womFileA.womType, Option(womFileA))
          .equals(
            WomOptionalValue(womFileB.womType, Option(womFileB))
          )
          .get
      }
    }

    it should s"fail symmetrically comparing (optionally) $description" in {
      a[WomExpressionException] should be thrownBy {
        WomOptionalValue(womFileB.womType, Option(womFileB))
          .equals(
            WomOptionalValue(womFileA.womType, Option(womFileA))
          )
          .get
      }
    }
  }

  val womStringEqualsTests = Table(
    ("description", "womFile", "string", "expected"),
    ("a single directory matched to a similar directory", singleDir, singleDir.value, true),
    ("a single file matched to a similar directory", singleFile, singleDir.value, false),
    ("a glob file matched to a similar directory", globFile, singleDir.value, false),
    ("a single directory matched to a similar file", singleDir, singleFile.value, false),
    ("a single file matched to a similar file", singleFile, singleFile.value, true),
    ("a glob file matched to a similar file", globFile, singleFile.value, false),
    ("a single directory matched to a similar glob", singleDir, globFile.value, false),
    ("a single file matched to a similar glob", singleFile, globFile.value, false),
    ("a glob file matched to a similar glob", globFile, globFile.value, true)
  )

  forAll(womStringEqualsTests) { (description, womFile, string, expected) =>
    it should s"expect $expected comparing $description as a string" in {
      womFile.equals(WomString(string)) should be(Success(WomBoolean(expected)))
    }

    it should s"expect $expected symmetrically comparing $description as a string" in {
      WomString(string).equals(womFile) should be(Success(WomBoolean(expected)))
    }

    it should s"expect $expected comparing $description as an optional string" in {
      WomOptionalValue(womFile.womType, Option(womFile)).equals(
        WomOptionalValue(WomStringType, Option(WomString(string)))
      ) should be(Success(WomBoolean(expected)))
    }

    it should s"expect $expected symmetrically comparing $description as an optional string" in {
      WomOptionalValue(WomStringType, Option(WomString(string))).equals(
        WomOptionalValue(womFile.womType, Option(womFile))
      ) should be(Success(WomBoolean(expected)))
    }
  }

  it should "produce an invalid equals" in {
    val result: Try[WomBoolean] = singleFile.equals(WomInteger(42))
    result.failed.get should have message "Cannot perform operation: single/file == WomInteger(42)"
  }
}
