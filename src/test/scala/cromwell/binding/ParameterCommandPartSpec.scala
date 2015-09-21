package cromwell.binding

import cromwell.binding.command.ParameterCommandPart
import cromwell.binding.types.{WdlIntegerType, WdlArrayType, WdlStringType}
import cromwell.binding.values.{WdlArray, WdlInteger, WdlString, WdlValue}
import cromwell.parser.BackendType
import cromwell.parser.WdlParser.SyntaxError
import org.scalatest.{FlatSpec, Matchers}

import scala.util.Failure

class ParameterCommandPartSpec extends FlatSpec with Matchers {
  val wdl =
    """task param_test {
      |  String a
      |  String b
      |  Array[String] c
      |  Int? d
      |  Array[Int] e
      |
      |  command <<<
      |  ./binary ${a} ${"-p " + b} ${sep="," c} ${default=9 d} ${sep="\t" e}
      |  >>>
      |}
      |
      |workflow wf {call param_test}
    """.stripMargin
  val namespace = WdlNamespace.load(wdl, BackendType.LOCAL)
  val task = namespace.tasks find {_.name == "param_test"} getOrElse {
    fail("task 'param_test' not found")
  }

  val paramsByName = task.commandTemplate.collect {case p: ParameterCommandPart => p}.map {p => p}
  "Template variables" should "Stringify correctly" in {
    paramsByName.size shouldEqual 5
    paramsByName(0).toString shouldEqual "${a}"
    paramsByName(1).toString shouldEqual "${\"-p \" + b}"
    paramsByName(2).toString shouldEqual "${sep=\",\" c}"
    paramsByName(3).toString shouldEqual "${default=\"9\" d}"
    paramsByName(4).toString shouldEqual "${sep=\"\\t\" e}"
  }

  "Command instantiation" should "succeed if given valid inputs" in {
    task.instantiateCommand(Map(
      "a" -> WdlString("a_val"),
      "b" -> WdlString("b_val"),
      "c" -> WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("c0"), WdlString("c1"), WdlString("c2"))),
      "d" -> WdlInteger(1),
      "e" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(0, 1, 2).map(WdlInteger(_)))
    )).get shouldEqual "./binary a_val -p b_val c0,c1,c2 1 0\t1\t2"
  }

  it should "succeed if omitting an optional input" in {
    task.instantiateCommand(Map(
      "a" -> WdlString("a_val"),
      "b" -> WdlString("b_val"),
      "c" -> WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("c0"), WdlString("c1"), WdlString("c2"))),
      "e" -> WdlArray(WdlArrayType(WdlIntegerType), Seq(0, 1, 2).map(WdlInteger(_)))
    )).get shouldEqual "./binary a_val -p b_val c0,c1,c2 9 0\t1\t2"
  }

  it should "succeed if providing an array with one element" in {
    task.instantiateCommand(Map(
      "a" -> WdlString("a_val"),
      "b" -> WdlString("b_val"),
      "c" -> WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("c0"))),
      "d" -> WdlInteger(1),
      "e" -> WdlArray(WdlArrayType(WdlIntegerType), Seq())
    )).get shouldEqual "./binary a_val -p b_val c0 1"
  }

  it should "succeed if providing an array with zero elements" in {
    task.instantiateCommand(Map(
      "a" -> WdlString("a_val"),
      "b" -> WdlString("b_val"),
      "c" -> WdlArray(WdlArrayType(WdlStringType), Seq()),
      "d" -> WdlInteger(1),
      "e" -> WdlArray(WdlArrayType(WdlIntegerType), Seq())
    )).get shouldEqual "./binary a_val -p b_val  1"
  }

  it should "raise exception if a required input is missing" in {
    task.instantiateCommand(Map("a" -> WdlString("a_val"))) match {
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
    )) match {
      case Failure(f) => // expected
      case _ => fail("Expected an exception")
    }
  }
}
