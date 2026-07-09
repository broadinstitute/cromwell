package wdl.shared.transforms.wdlom2wom

import cats.data.Validated.{Invalid, Valid}
import org.scalatest.flatspec.AnyFlatSpec
import wom.types.{WomIntegerType, WomObjectType, WomStringType, WomType}
import wom.values.{WomInteger, WomObject, WomString, WomValue}

class WdlSharedInputParsingSpec extends AnyFlatSpec {

  behavior of "WdlSharedInputParsing"

  it should "parse inputs without runtime overrides correctly" in {
    val inputJson =
      """
        |{
        |  "wf.task1.input1": "value1",
        |  "wf.task2.input2": 42
        |}
        |""".stripMargin

    val expectedParsedInputs = Map(
      "wf.task1.input1" -> WomString("value1"),
      "wf.task2.input2" -> WomInteger(42)
    )

    val typesMap = Map(
      "wf.task1.input1" -> WomStringType,
      "wf.task2.input2" -> WomIntegerType
    )

    check(inputJson, expectedParsedInputs, typesMap)

  }

  it should "parse inputs with a single runtime override correctly" in {
    val inputJson =
      """
        |{
        |  "wf.task1.input1": "value1",
        |  "wf.task2.input2": 42,
        |  "wf.task1.runtime.cpu": "2"
        |}
        |""".stripMargin

    val expectedParsedInputs = Map(
      "wf.task1.input1" -> WomString("value1"),
      "wf.task2.input2" -> WomInteger(42),
      "wf.task1.runtime" -> WomObject(Map("cpu" -> WomString("2")))
    )

    val typeMap = Map(
      "wf.task1.input1" -> WomStringType,
      "wf.task2.input2" -> WomIntegerType,
      "wf.task1.runtime" -> WomObjectType
    )

    check(inputJson, expectedParsedInputs, typeMap)

  }

  it should "parse inputs with multiple runtime overrides correctly" in {
    val inputJson =
      """
        |{
        |  "wf.task1.input1": "value1",
        |  "wf.task2.input2": 42,
        |  "wf.task1.runtime.cpu": "2",
        |  "wf.task1.runtime.memory": "4 GB",
        |  "wf.task1.runtime.disk": "10 GB",
        |  "wf.task2.runtime.cpu": "4",
        |  "wf.task2.runtime.memory": "8 GB",
        |  "wf.task2.input2": 42,
        |  "wf.sub_wf.task3.runtime.gpu": "foo",
        |  "wf.sub_wf.task3.input3": "bar"
        |}
        |""".stripMargin

    val expectedParsedInputs = Map(
      "wf.task1.input1" -> WomString("value1"),
      "wf.task1.runtime" -> WomObject(
        Map(
          "cpu" -> WomString("2"),
          "memory" -> WomString("4 GB"),
          "disk" -> WomString("10 GB")
        )
      ),
      "wf.task2.input2" -> WomInteger(42),
      "wf.task2.runtime" -> WomObject(
        Map(
          "cpu" -> WomString("4"),
          "memory" -> WomString("8 GB")
        )
      ),
      "wf.sub_wf.task3.input3" -> WomString("bar"),
      "wf.sub_wf.task3.runtime" -> WomObject(
        Map(
          "gpu" -> WomString("foo")
        )
      )
    )

    val typeMap = Map(
      "wf.task1.input1" -> WomStringType,
      "wf.task2.input2" -> WomIntegerType,
      "wf.sub_wf.task3.input3" -> WomStringType,
      "wf.task1.runtime" -> WomObjectType,
      "wf.task2.runtime" -> WomObjectType,
      "wf.sub_wf.task3.runtime" -> WomObjectType
    )

    check(inputJson, expectedParsedInputs, typeMap)

  }

  def check(jsonString: String, expected: Map[String, WomValue], inputsToTypes: Map[String, WomType]) =
    WdlSharedInputParsing.parseInputs(jsonString) match {
      case Right(parsedInputs) =>
        val coercedInputs = parsedInputs.map { case (k, v) =>
          val coercedInput = v(inputsToTypes(k)) match {
            case Valid(coerced) => coerced
            case Invalid(error) => fail(s"Coercion failed for key $k with error: $error")
          }
          k -> coercedInput
        }
        assert(coercedInputs == expected)
      case Left(error) => fail(s"Parsing failed with error: $error")
    }
}
