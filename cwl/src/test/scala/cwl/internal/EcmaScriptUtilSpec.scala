package cwl.internal

import common.validation.Validation._
import delight.rhinosandox.exceptions.ScriptDurationException
import org.mozilla.javascript.{EcmaError, EvaluatorException}
import org.scalatest.{FlatSpec, Matchers}
import wom.types._
import wom.values._

class EcmaScriptUtilSpec extends FlatSpec with Matchers {

  behavior of "EcmaScriptUtil"

  private lazy val encoder = new EcmaScriptEncoder()

  it should "lookup a map entry with a string key" in {
    val values =
      "myName" -> WomMap(
        WomMapType(WomBooleanType, WomArrayType(WomStringType)),
        Map(WomBoolean(true) -> WomArray(WomArrayType(WomStringType), Seq(WomString("myValue"))))
      )
    val expr = """myName["true"][0] + 'Plus'"""
    val result: WomValue = EcmaScriptUtil.evalStructish(expr, values, encoder = encoder).toTry.get
    result should be(WomString("myValuePlus"))
  }

  it should "lookup a map entry with a boolean key" in {
    val values =
      "myName" -> WomMap(
        WomMapType(WomBooleanType, WomArrayType(WomStringType)),
        Map(WomBoolean(true) -> WomArray(WomArrayType(WomStringType), Seq(WomString("myValue"))))
      )
    val expr = "myName[true][0] + 'PlusPlus'"
    val result = EcmaScriptUtil.evalStructish(expr, values, encoder = encoder).toTry.get
    result should be(WomString("myValuePlusPlus"))
  }

  it should "JSON.stringify" in {
    val values =
      "myName" -> WomMap(
        WomMapType(WomBooleanType, WomArrayType(WomStringType)),
        Map(WomBoolean(true) -> WomArray(WomArrayType(WomStringType), Seq(WomString("myValue"))))
      )
    val expr = "JSON.stringify(myName)"
    val result: WomValue = EcmaScriptUtil.evalStructish(expr, values, encoder = encoder).toTry.get
    result should be(WomString("""{"true":["myValue"]}"""))
  }

  it should "Array.prototype.sort" in {
    val values = "myName" -> WomArray(List("hello", "world", "alpha", "zulu", "indigo").map(WomString))
    val expr = "myName.sort()"
    val expected = WomArray(List("alpha", "hello", "indigo", "world", "zulu").map(WomString))
    val result = EcmaScriptUtil.evalStructish(expr, values, encoder = encoder).toTry.get
    result should be(expected)
  }

  it should "run in strict mode" in {
    the[EvaluatorException] thrownBy {
      EcmaScriptUtil.evalRaw(
        """|function sum(a, a, c) { // !!! syntax error
           |  return a + a + c; // wrong if this code ran
           |}
           |""".stripMargin) { (_, _) => () }
    } should have message
      """Parameter "a" already declared in this function. (<ecmascript>#1)"""
  }

  it should "not run java code" in {
    the[EcmaError] thrownBy {
      EcmaScriptUtil.evalRaw(
        """|var path = java.nio.file.Paths.get("/etc/passwd");
           |java.nio.file.Files.exists(path);
           |// Or even... java.lang.System.exit(42);
           |""".stripMargin) { (_, _) => () }
    } should have message
      """ReferenceError: "java" is not defined. (<ecmascript>#1)"""
  }

  it should "not run forever" in {
    a[ScriptDurationException] should be thrownBy {
      EcmaScriptUtil.evalRaw(
        """|function sleep(milliseconds) {
           |  var start = new Date().getTime();
           |  while(true) {
           |    if ((new Date().getTime() - start) > milliseconds){
           |      break;
           |    }
           |  }
           |}
           |sleep(120 * 1000);
           |""".stripMargin) { (_, _) => () }
    }
  }
}
