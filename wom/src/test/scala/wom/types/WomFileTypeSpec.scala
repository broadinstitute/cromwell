package wom.types

import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import spray.json.{JsNumber, JsString}
import wom.WomExpressionException
import wom.values.{WomFloat, WomGlobFile, WomSingleFile, WomString, WomUnlistedDirectory}

import scala.util.Success

class WomFileTypeSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "WomFileType"

  lazy val coercionTests = Table(
    ("description", "value", "expected"),

    ("a string to a dir", "example/string", WomUnlistedDirectory("example/string")),
    ("a string to a file", "example/string", WomSingleFile("example/string")),
    ("a string to a glob", "example/string", WomGlobFile("example/string")),

    ("a js string to a dir", JsString("example/js"), WomUnlistedDirectory("example/js")),
    ("a js string to a file", JsString("example/js"), WomSingleFile("example/js")),
    ("a js string to a glob", JsString("example/js"), WomGlobFile("example/js")),

    ("a wom string to a dir", WomString("example/wom"), WomUnlistedDirectory("example/wom")),
    ("a wom string to a file", WomString("example/wom"), WomSingleFile("example/wom")),
    ("a wom string to a glob", WomString("example/wom"), WomGlobFile("example/wom")),

    ("a wom dir to a dir", WomUnlistedDirectory("example/dir"),
      WomUnlistedDirectory("example/dir")),
    ("a wom file to a file", WomSingleFile("example/dir"), WomSingleFile("example/dir")),
    ("a wom glob to a glob", WomGlobFile("example/glob*"), WomGlobFile("example/glob*")),

    ("a string with leading and trailing special symbols to a dir", "  \n\r\t example/string \n\r\t ",
      WomUnlistedDirectory("example/string")),
    ("a string with leading and trailing special symbols to a file", "  \n\r\t example/string \n\r\t ",
      WomSingleFile("example/string")),
    ("a string with leading and trailing special symbols to a glob", "  \n\r\t example/string \n\r\t ",
      WomGlobFile("example/string")),

    ("a js string with leading and trailing special symbols to a dir", JsString("  \n\r\t example/js \n\r\t "),
      WomUnlistedDirectory("example/js")),
    ("a js string with leading and trailing special symbols to a file", JsString("  \n\r\t example/js \n\r\t "),
      WomSingleFile("example/js")),
    ("a js string with leading and trailing special symbols to a glob", JsString("  \n\r\t example/js \n\r\t "),
      WomGlobFile("example/js")),

    ("a wom string with leading and trailing special symbols to a dir", WomString("  \n\r\t example/wom \n\r\t "),
      WomUnlistedDirectory("example/wom")),
    ("a wom string with leading and trailing special symbols to a file", WomString("  \n\r\t example/wom \n\r\t "),
      WomSingleFile("example/wom")),
    ("a wom string with leading and trailing special symbols to a glob", WomString("  \n\r\t example/wom \n\r\t "),
      WomGlobFile("example/wom")),
  )

  forAll(coercionTests) { (description, value, expected) =>
    it should s"coerce $description" in {
      val womFileType = expected.womType
      womFileType.coerceRawValue(value).get should be(expected)
      womFileType.coercionDefined(value) should be(true)
    }
  }

  lazy val failedCoercionTests = Table(
    ("description", "womFileType", "value", "expected"),

    ("a double to a dir", WomUnlistedDirectoryType, 6.28318,
      "No coercion defined from '6.28318' of type 'java.lang.Double' to 'Directory'."),
    ("a double to a file", WomSingleFileType, 6.28318,
      "No coercion defined from '6.28318' of type 'java.lang.Double' to 'File'."),
    ("a double to a glob", WomGlobFileType, 6.28318,
      "No coercion defined from '6.28318' of type 'java.lang.Double' to 'Glob'."),

    ("a js number to a dir", WomUnlistedDirectoryType, JsNumber(6.28318),
      "No coercion defined from '6.28318' of type 'spray.json.JsNumber' to 'Directory'."),
    ("a js number to a file", WomSingleFileType, JsNumber(6.28318),
      "No coercion defined from '6.28318' of type 'spray.json.JsNumber' to 'File'."),
    ("a js number to a glob", WomGlobFileType, JsNumber(6.28318),
      "No coercion defined from '6.28318' of type 'spray.json.JsNumber' to 'Glob'."),

    ("a wom float to a dir", WomUnlistedDirectoryType, WomFloat(6.28318),
      "No coercion defined from wom value(s) '6.28318' of type 'Float' to 'Directory'."),
    ("a wom float to a file", WomSingleFileType, WomFloat(6.28318),
      "No coercion defined from wom value(s) '6.28318' of type 'Float' to 'File'."),
    ("a wom float to a glob", WomGlobFileType, WomFloat(6.28318),
      "No coercion defined from wom value(s) '6.28318' of type 'Float' to 'Glob'."),

    ("a wom file to a dir", WomUnlistedDirectoryType, WomSingleFile("example/file"),
      """No coercion defined from wom value(s) '"example/file"' of type 'File' to 'Directory'."""),
    ("a wom glob to a dir", WomUnlistedDirectoryType, WomGlobFile("example/glob*"),
      """No coercion defined from wom value(s) 'glob("example/glob*")' of type 'Glob' to 'Directory'."""),

    ("a wom dir to a file", WomSingleFileType, WomUnlistedDirectory("example/dir"),
      """No coercion defined from wom value(s) '"example/dir"' of type 'Directory' to 'File'."""),
    ("a wom glob to a file", WomSingleFileType, WomGlobFile("example/glob*"),
      """No coercion defined from wom value(s) 'glob("example/glob*")' of type 'Glob' to 'File'."""),

    ("a wom dir to a glob", WomGlobFileType, WomUnlistedDirectory("example/dir"),
      """No coercion defined from wom value(s) '"example/dir"' of type 'Directory' to 'Glob'."""),
    ("a wom file to a glob", WomGlobFileType, WomSingleFile("example/file"),
      """No coercion defined from wom value(s) '"example/file"' of type 'File' to 'Glob'.""")
  )

  forAll(failedCoercionTests) { (description, womFileType, value, expected) =>
    it should s"fail to coerce $description" in {
      the[IllegalArgumentException] thrownBy {
        womFileType.coerceRawValue(value).get
      } should have message expected
      womFileType.coercionDefined(value) should be(false)
    }
  }

  lazy val womFileTypes = Table(
    ("womFileTypeName", "womFileType"),
    ("WomUnlistedDirectoryType", WomUnlistedDirectoryType),
    ("WomSingleFileType", WomSingleFileType),
    ("WomGlobFileType", WomGlobFileType)
  )

  forAll(womFileTypes) { (womFileTypeName, womFileType) =>
    it should s"add a $womFileTypeName with a WomStringType" in {
      womFileType.add(WomStringType) should be(Success(womFileType))
    }

    it should s"add a WomStringType with a $womFileTypeName" in {
      WomStringType.add(womFileType) should be(Success(WomStringType))
    }

    it should s"add an optional $womFileTypeName with an optional WomStringType" in {
      WomOptionalType(womFileType).add(WomOptionalType(WomStringType)) should be(Success(womFileType))
    }

    it should s"add an optional WomStringType with an optional $womFileTypeName" in {
      WomOptionalType(WomStringType).add(WomOptionalType(womFileType)) should be(Success(WomStringType))
    }

    it should s"not add a $womFileTypeName with a WomFloatType" in {
      the[WomExpressionException] thrownBy {
        womFileType.add(WomFloatType).get
      } should have message s"Type evaluation cannot determine type from expression: $womFileTypeName + WomFloatType"
    }

    it should s"not add an optional $womFileTypeName with a optional WomFloatType" in {
      the[WomExpressionException] thrownBy {
        WomOptionalType(womFileType).add(WomOptionalType(WomFloatType)).get
      } should have message s"Type evaluation cannot determine type from expression: $womFileTypeName + WomFloatType"
    }

    it should s"equal a $womFileTypeName with a WomStringType" in {
      womFileType.equalsType(WomStringType) should be(Success(WomBooleanType))
    }

    it should s"equal an optional $womFileTypeName with an optional WomStringType" in {
      WomOptionalType(womFileType).equalsType(WomOptionalType(WomStringType)) should be(Success(WomBooleanType))
    }

    it should s"not compare a $womFileTypeName with a WomFloatType" in {
      the[WomExpressionException] thrownBy {
        womFileType.equalsType(WomFloatType).get
      } should have message s"Type evaluation cannot determine type from expression: $womFileTypeName == WomFloatType"
    }

    it should s"not compare an optional $womFileTypeName with an optional WomFloatType" in {
      the[WomExpressionException] thrownBy {
        WomOptionalType(womFileType).equalsType(WomOptionalType(WomFloatType)).get
      } should have message s"Type evaluation cannot determine type from expression: $womFileTypeName == WomFloatType"
    }
  }
}
