package cwl

import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{BeforeAndAfterEach, FlatSpec, Matchers}
import org.specs2.mock.Mockito
import spray.json.{DefaultJsonProtocol, JsArray, JsBoolean, JsNull, JsNumber, JsObject, JsString, JsValue, RootJsonFormat}
import wom.graph.ScatterNode.ScatterVariableAndValue
import wom.graph.expression.PlainAnonymousExpressionNode
import wom.graph.{ScatterNode, ScatterVariableNode, WomIdentifier}
import wom.types.{WomArrayType, WomStringType}
import wom.values.{WomArray, WomBoolean, WomFile, WomFloat, WomInteger, WomMap, WomObject, WomOptionalValue, WomPair, WomString, WomValue}

class ScatterLogicSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks with Mockito  with DefaultJsonProtocol with BeforeAndAfterEach {
  private val expressionNode = PlainAnonymousExpressionNode(WomIdentifier("name"), null, WomStringType, Map.empty)
  private val arrayStringType = WomArrayType(WomStringType)

  val emptyArray = ScatterVariableAndValue(
    ScatterVariableNode(null, expressionNode, WomStringType),
    WomArray(arrayStringType, List.empty)
  )
  val simpleArray2 = ScatterVariableAndValue(
    ScatterVariableNode(null, expressionNode, WomStringType),
    WomArray(arrayStringType, List(WomString("a"), WomString("b")))
  )
  val simpleArray3 = ScatterVariableAndValue(
    ScatterVariableNode(null, expressionNode, WomStringType),
    WomArray(arrayStringType, List(WomString("a"), WomString("b"), WomString("c")))
  )
  val simpleArray4 = ScatterVariableAndValue(
    ScatterVariableNode(null, expressionNode, WomStringType),
    WomArray(arrayStringType, List(WomString("a"), WomString("b"), WomString("c"), WomString("d")))
  )
  
  override def beforeEach() = {
    // The index length of the SVN is a mutable value, to avoid tests stepping on each others reset it to the default value
    // before each test
    emptyArray.scatterVariableNode.withRelativeIndexLength(1)
    simpleArray2.scatterVariableNode.withRelativeIndexLength(1)
    simpleArray3.scatterVariableNode.withRelativeIndexLength(1)
    simpleArray4.scatterVariableNode.withRelativeIndexLength(1)
  }

  val valueA = WomString("a")
  val valueB = WomString("b")
  val valueC = WomString("c")
  val valueD = WomString("d")
  
  def validateScatterCombinations(list: List[ScatterNode.ScatterVariableAndValue],
                                  scatterSize: Int,
                                  expected: List[List[Int]]) = {
    // Go through all the variable nodes and feed all the shard numbers to the indexForShard function
    val combinations = list.map({
      case ScatterVariableAndValue(node, values) =>
        (0 until scatterSize).map(node.indexForShard(_, values.size)).toList
    })

    combinations should contain theSameElementsInOrderAs expected
  }

  "CrossProduct scattering function" should "generate correct results" in {
    val nonEmptyVariables = Table(
      // None for the index length means "don't check because irrelevant"
      ("variablesAndValues", "scatterSize", "arrayIndexes"),
      (List(simpleArray2), 2, List(List(0, 1))),
      (List(simpleArray2, simpleArray3), 6, List(List(0, 0, 0, 1, 1, 1), List(0, 1, 2, 0, 1, 2))),
      (List(simpleArray2, simpleArray3, simpleArray4), 24, List(
        List(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
        List(0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2),
        List(0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3))
      )
    )

    forAll(nonEmptyVariables) { case (list, scatterSize, expected) =>
      ScatterLogic.CrossProductScatterProcessingFunction(list) shouldBe Right(scatterSize)
      validateScatterCombinations(list, scatterSize, expected)
    }

    List(List(emptyArray), List(emptyArray, simpleArray2), List(simpleArray2, emptyArray)) foreach { list =>
      ScatterLogic.CrossProductScatterProcessingFunction(list) shouldBe Right(0)
    }
  }

  "DotProduct scattering function" should "generate correct results" in {
    import common.validation.Checked._
    val simpleArray2Duplicate = ScatterVariableAndValue(
      ScatterVariableNode(null, expressionNode, WomStringType),
      WomArray(arrayStringType, List(WomString("a"), WomString("b")))
    )
    val nonEmptyVariables = Table(
      // None for the index length means "don't check because irrelevant"
      ("variablesAndValues", "scatterSize", "arrayIndexes"),
      (List(simpleArray2, simpleArray2Duplicate), 2.validNelCheck, List(List(0, 1), List(0, 1))),
      (List(simpleArray2, simpleArray3), "All arrays must have the same number of element when using the dot product scatter method".invalidNelCheck, List.empty)
    )

    forAll(nonEmptyVariables) { case (list, expectedScatterSize, expected) =>
      ScatterLogic.DotProductScatterProcessingFunction(list) shouldBe expectedScatterSize

      expectedScatterSize match {
        case Right(scatterSize) => validateScatterCombinations(list, scatterSize, expected)
        case _ => // Only verify combinations if the scatter is valid
      }
    }
  }

  "NestedCrossProduct collecting function builder" should "generate correct results" in {
    import spray.json._

    val arrayString = arrayStringType
    val nestedArrayString = WomArrayType(arrayString)
    val nestedArrayString2 = WomArrayType(nestedArrayString)

    val variables = Table(
      ("arraySizes", "shards", "valueType", "result"),
      (List(0), List.empty, arrayString, "[]"),
      (List(0, 2), List.empty, nestedArrayString, "[]"),
      (List(2, 0), List.empty, nestedArrayString, "[[], []]"),
      (List(2, 3, 0), List.empty, nestedArrayString2, "[[[], [], []], [[], [], []]]"),
      (List(2, 0, 3), List.empty, nestedArrayString2, "[[], []]"),
      // This technically should not happen since to have a cross product we need at least 2 scatter variables but..
      (List(2), List(valueA, valueB), arrayString, """["a", "b"]"""),
      (List(2, 2), List(valueA, valueB, valueC, valueD), nestedArrayString, """[["a", "b"], ["c", "d"]]"""),
      (List(3, 2), List(valueA, valueB, valueC, valueD, valueA, valueB), nestedArrayString, """[["a", "b"], ["c", "d"], ["a", "b"]]"""),
      (List(2, 3), List(valueA, valueB, valueC, valueC, valueB, valueA), nestedArrayString, """[["a", "b", "c"], ["c", "b", "a"]]"""),
      (List(2, 3, 2), List(
        valueA, valueB, valueA, valueB, valueA, valueB,
        valueC, valueD, valueC, valueD, valueC, valueD
      ), nestedArrayString2, """[[["a", "b"], ["a", "b"], ["a", "b"]], [["c", "d"], ["c", "d"], ["c", "d"]]]""")
    )

    forAll(variables) { case (arraySizes, shards, valueType, result) =>
      val collection = ScatterLogic.NestedCrossProductScatterCollectionFunctionBuilder(arraySizes)(shards, valueType)
      collection.arrayType shouldBe valueType
      (collection: WomValue).toJson(WomValueJsonFormatter.WomValueJsonFormat) shouldBe result.parseJson
    }
  }

  // This exist in engine currently but can't be used from here..
  object WomValueJsonFormatter {
    implicit object WomValueJsonFormat extends RootJsonFormat[WomValue] {
      def write(value: WomValue): JsValue = value match {
        case s: WomString => JsString(s.value)
        case i: WomInteger => JsNumber(i.value)
        case f: WomFloat => JsNumber(f.value)
        case b: WomBoolean => JsBoolean(b.value)
        case f: WomFile => JsString(f.value)
        case o: WomObject => new JsObject(o.value map {case(k, v) => k -> write(v)})
        case a: WomArray => new JsArray(a.value.map(write).toVector)
        case m: WomMap => new JsObject(m.value map {case(k,v) => k.valueString -> write(v)})
        case q: WomPair => new JsObject(Map("left" -> write(q.left), "right" -> write(q.right)))
        case WomOptionalValue(_, Some(innerValue)) => write(innerValue)
        case WomOptionalValue(_, None) => JsNull
        // handles WdlExpression
        case v: WomValue => JsString(v.toWomString)

      }
      def read(value: JsValue): WomValue = ???
    }
  }
}
