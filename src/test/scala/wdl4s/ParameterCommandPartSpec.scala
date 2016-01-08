package wdl4s

import wdl4s.command.ParameterCommandPart
import wdl4s.expression.NoFunctions
import wdl4s.types.{WdlIntegerType, WdlArrayType, WdlStringType}
import wdl4s.values._
import wdl4s.parser.WdlParser.SyntaxError
import org.scalatest.{FlatSpec, Matchers}

import scala.util.{Try, Failure}

class ParameterCommandPartSpec extends FlatSpec with Matchers {
  val wdl =
    """task param_test {
      |  String a
      |  String b
      |  Array[String] c
      |  Int? d
      |  Array[Int] e
      |  Boolean f
      |
      |  command <<<
      |  ./binary ${a} ${"-p " + b} ${sep="," c} ${default=9 d} ${sep="\t" e} ${true="--true" false="--false" f}
      |  >>>
      |}
      |
      |workflow wf {call param_test}
    """.stripMargin
  val namespace = WdlNamespace.load(wdl)
  val task = namespace.tasks find {_.name == "param_test"} getOrElse {
    fail("task 'param_test' not found")
  }

  val paramsByName = task.commandTemplate.collect {case p: ParameterCommandPart => p}.map {p => p}
  "Template variables" should "Stringify correctly" in {
    paramsByName.size shouldEqual 6
    paramsByName.head.toString shouldEqual "${a}"
    paramsByName(1).toString shouldEqual "${\"-p \" + b}"
    paramsByName(2).toString shouldEqual "${sep=\",\" c}"
    paramsByName(3).toString shouldEqual "${default=\"9\" d}"
    paramsByName(4).toString shouldEqual "${sep=\"\\t\" e}"
    paramsByName(5).toString shouldEqual "${true=\"--true\" false=\"--false\" f}"
  }

  "Command instantiation" should "succeed if given valid inputs" in {
    task.instantiateCommand(Map(
      "a" -> WdlString("a_val"),
      "b" -> WdlString("b_val"),
      "c" -> WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("c0"), WdlString("c1"), WdlString("c2"))),
      "d" -> WdlInteger(1),
      "e" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(0, 1, 2).map(WdlInteger(_))),
      "f" -> WdlBoolean.False
    ), NoFunctions).get shouldEqual "./binary a_val -p b_val c0,c1,c2 1 0\t1\t2 --false"
  }

  it should "succeed if omitting an optional input" in {
    task.instantiateCommand(Map(
      "a" -> WdlString("a_val"),
      "b" -> WdlString("b_val"),
      "c" -> WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("c0"), WdlString("c1"), WdlString("c2"))),
      "e" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(0, 1, 2).map(WdlInteger(_))),
      "f" -> WdlBoolean.True
    ), NoFunctions).get shouldEqual "./binary a_val -p b_val c0,c1,c2 9 0\t1\t2 --true"
  }

  it should "succeed if providing an array with one element" in {
    task.instantiateCommand(Map(
      "a" -> WdlString("a_val"),
      "b" -> WdlString("b_val"),
      "c" -> WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("c0"))),
      "d" -> WdlInteger(1),
      "e" -> WdlArray(WdlArrayType(WdlIntegerType), Seq()),
      "f" -> WdlBoolean.True
    ), NoFunctions).get shouldEqual "./binary a_val -p b_val c0 1  --true"
  }

  it should "succeed if providing an array with zero elements" in {
    task.instantiateCommand(Map(
      "a" -> WdlString("a_val"),
      "b" -> WdlString("b_val"),
      "c" -> WdlArray(WdlArrayType(WdlStringType), Seq()),
      "d" -> WdlInteger(1),
      "e" -> WdlArray(WdlArrayType(WdlIntegerType), Seq()),
      "f" -> WdlBoolean.True
    ), NoFunctions).get shouldEqual "./binary a_val -p b_val  1  --true"
  }

  it should "raise exception if a required input is missing" in {
    task.instantiateCommand(Map("a" -> WdlString("a_val")), NoFunctions) match {
      case Failure(f) => // expected
      case _ => fail("Expected an exception")
    }
  }

  it should "raise exception if a parameter is an expression" in {
    task.instantiateCommand(Map(
      "a" -> WdlString("a_val"),
      "b" -> WdlExpression.fromString("'a'+'b'"),
      "c" -> WdlArray(WdlArrayType(WdlStringType), Seq()),
      "d" -> WdlInteger(1),
      "e" -> WdlArray(WdlArrayType(WdlIntegerType), Seq())
    ), NoFunctions) match {
      case Failure(f) => // expected
      case _ => fail("Expected an exception")
    }
  }

  it should "raise exception if 'true' attribute is specified but 'false' is not" in {
    Try(
      WdlNamespace.load(
        """task param_test {
          |  Boolean f
          |
          |  command <<<
          |  ./binary ${true="--true" f}
          |  >>>
          |}
          |
          |workflow wf {call param_test}
        """.stripMargin)
    ) match {
      case Failure(s: SyntaxError) => // expected
      case _ => fail("Expecting a syntax error")
    }
  }
}
