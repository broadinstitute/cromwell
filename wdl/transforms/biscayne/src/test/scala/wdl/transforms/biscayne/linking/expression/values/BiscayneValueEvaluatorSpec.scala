package wdl.transforms.biscayne.linking.expression.values

import common.assertion.CromwellTimeoutSpec
import common.assertion.ErrorOrAssertions._
import common.assertion.ManyTimes.intWithTimes
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.expression.EvaluatedValue
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wdl.transforms.biscayne.Ast2WdlomSpec.{fromString, parser}
import wdl.transforms.biscayne.ast2wdlom._
import wom.expression.NoIoFunctionSet
import wom.types.{WomAnyType, WomArrayType, WomIntegerType, WomMapType, WomOptionalType, WomStringType}
import wom.values.{WomArray, WomBoolean, WomInteger, WomMap, WomObject, WomOptionalValue, WomPair, WomString}

class BiscayneValueEvaluatorSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  behavior of "biscayne value evaluator"

  // The "did we bring in the base transforms" sanity check:
  it should "return a static integer from static integer addition" in {
    val str = "3 + 3"
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(WomInteger(6), Seq.empty)
    }
  }

  it should "evaluate an 'as_map' expression correctly" in {
    val str = """ as_map( [("x", 1), ("y", 2), ("z", 3)] ) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedMap: WomMap = WomMap(
      Map(
        WomString("x") -> WomInteger(1),
        WomString("y") -> WomInteger(2),
        WomString("z") -> WomInteger(3)
      )
    )

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedMap, Seq.empty)
    }
  }

  it should "evaluate an 'as_pairs' expression correctly" in {
    // Repeat this a few times ensure we aren't occasionally reordering by accident:
    10 times {
      val str = """ as_pairs( { 1: "one", 2: "two", 3: three } ) """
      val expr = fromString[ExpressionElement](str, parser.parse_e)

      val inputs = Map("three" -> WomString("three"))
      val expectedPairs: WomArray = WomArray(
        Seq(
          WomPair(WomInteger(1), WomString("one")),
          WomPair(WomInteger(2), WomString("two")),
          WomPair(WomInteger(3), WomString("three"))
        )
      )

      expr.shouldBeValidPF { case e =>
        e.evaluateValue(inputs, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedPairs, Seq.empty)
      }
      ()
    }
  }

  // TODO: Sort out stable map ordering
  // We cannot run this test yet because the map reorders the pairs
  it should "echo correctly via as_pairs(as_map(...))" ignore {

    // Repeat this a few times ensure we aren't occasionally reordering by accident:
    10 times {
      // A value lookup expression:
      val str = """ as_pairs(as_map(echo_me)) """
      val expr = fromString[ExpressionElement](str, parser.parse_e)

      val expectedPairs: WomArray = WomArray(
        Seq(
          WomPair(WomInteger(1), WomString("one")),
          WomPair(WomInteger(2), WomString("two")),
          WomPair(WomInteger(3), WomString("three"))
        )
      )

      val inputs = Map("echo_me" -> expectedPairs)

      expr.shouldBeValidPF { case e =>
        e.evaluateValue(inputs, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedPairs, Seq.empty)
      }
      ()
    }
  }

  it should "fail to evaluate 'as_map' if keys are duplicated" in {

    // Repeat this a few times ensure we aren't occasionally reordering by accident:
    10 times {
      val str = """ as_map( [("x", 1), ("y", 2), ("x", 3)] ) """
      val expr = fromString[ExpressionElement](str, parser.parse_e)

      expr.shouldBeValidPF { case e =>
        e.evaluateValue(Map.empty, NoIoFunctionSet, None)
          .shouldBeInvalid(
            """Cannot evaluate 'as_map' with duplicated keys: keys can only appear once but "x" appeared 2 times."""
          )
      }
      ()
    }
  }

  it should "collect by key when keys are duplicated" in {

    // Repeat this a few times ensure we aren't occasionally reordering by accident:
    10 times {
      val str = """ collect_by_key( [("x", 1), ("y", 2), ("x", 3)] ) """
      val expr = fromString[ExpressionElement](str, parser.parse_e)

      val expectedMap: WomMap = WomMap(
        Map(
          WomString("x") -> WomArray(Seq(WomInteger(1), WomInteger(3))),
          WomString("y") -> WomArray(Seq(WomInteger(2)))
        )
      )

      expr.shouldBeValidPF { case e =>
        e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedMap, Seq.empty)
      }
      ()
    }
  }

  it should "evaluate a map literal mixing String?s and Int?s" in {
    val str = """ { "i": i_in, "s": s_in } """
    val inputs = Map(
      "i_in" -> WomOptionalValue(WomIntegerType, Some(WomInteger(1))),
      "s_in" -> WomOptionalValue(WomStringType, Some(WomString("two")))
    )
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedMap: WomMap = WomMap(
      WomMapType(WomStringType, WomOptionalType(WomStringType)),
      Map(
        WomString("i") -> WomOptionalValue(WomStringType, Some(WomString("1"))),
        WomString("s") -> WomOptionalValue(WomStringType, Some(WomString("two")))
      )
    )

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(inputs, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedMap, Seq.empty)
    }
  }

  val escapeTests = Map(
    "\\\\" -> "\\",
    "\\n" -> System.lineSeparator,
    "\\t" -> "\t",
    "\\'" -> "'",
    "\\\"" -> "\"",
    "\\150\\145\\154\\154\\157" -> "hello",
    "\\x68\\x65\\x6C\\x6c\\x6F" -> "hello",
    "\\u0068\\U00000065\\u006C\\U0000006C\\u006F" -> "hello",
    "\\u03A9 (omega)" -> "Ω (omega)"
  )

  escapeTests foreach { case (sequence, expected) =>
    List("\"", "'") foreach { quote =>
      it should s"evaluate the escaping string $quote$sequence$quote as: $expected" in {
        val str = s"$quote$sequence$quote"
        val expectedEvaluation = WomString(expected)
        val expr = fromString[ExpressionElement](str, parser.parse_e)

        expr.shouldBeValidPF { case e =>
          e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedEvaluation, Seq.empty)
        }
      }
    }
  }

  it should "evaluate a simple sep expression correctly" in {
    val str = """ sep(" ", ["a", "b", "c"]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedString: WomString = WomString("a b c")

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedString, Seq.empty)
    }
  }

  it should "evaluate a sep expression containing a sub-call to prefix correctly" in {
    val str = """ sep(" ", prefix("-i ", ["a", "b", "c"])) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedString: WomString = WomString("-i a -i b -i c")

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedString, Seq.empty)
    }
  }

  it should "evaluate a POSIX-flavor regex in a sub expression correctly" in {
    val str = """ sub("aB", "[[:lower:]]", "9") """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedString: WomString = WomString("9B")

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedString, Seq.empty)
    }
  }

  it should "fail to evaluate a Java-flavor regex in a sub expression" in {
    val str = """ sub("aB", "\\p{Lower}", "9") """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None)
        .shouldBeInvalid(
          """error parsing regexp: invalid character class range: `\p{Lower}`"""
        )
    }
  }

  it should "evaluate a suffix expression correctly" in {
    val str = """ suffix("S", ["a", "b", "c"]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedArray: WomArray = WomArray(
      Seq(
        WomString("aS"),
        WomString("bS"),
        WomString("cS")
      )
    )

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedArray, Seq.empty)
    }
  }

  it should "evaluate a quote expression correctly with an empty array" in {
    val str = """ quote([]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedArray: WomArray = WomArray(Seq())

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedArray, Seq.empty)
    }
  }

  it should "evaluate a quote expression correctly with an array of integers" in {
    val str = """ quote([1, 2, 3]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedArray: WomArray = WomArray(
      Seq(
        WomString("\"1\""),
        WomString("\"2\""),
        WomString("\"3\"")
      )
    )

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedArray, Seq.empty)
    }
  }

  it should "evaluate a quote expression correctly with an array of booleans" in {
    val str = """ quote([true, false, true]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedArray: WomArray = WomArray(
      Seq(
        WomString("\"true\""),
        WomString("\"false\""),
        WomString("\"true\"")
      )
    )

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedArray, Seq.empty)
    }
  }

  it should "evaluate a quote expression correctly with an array of floats" in {
    val str = """ quote([1.1, 2.2, 3]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedArray: WomArray = WomArray(
      Seq(
        WomString("\"1.1\""),
        WomString("\"2.2\""),
        WomString("\"3.0\"")
      )
    )

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedArray, Seq.empty)
    }
  }

  it should "evaluate a quote expression correctly with an array of files" in {
    val str = """ quote(["/home/someFile.txt", "/rootFile.exe", "./anotherFile.py"]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedArray: WomArray = WomArray(
      Seq(
        WomString("\"/home/someFile.txt\""),
        WomString("\"/rootFile.exe\""),
        WomString("\"./anotherFile.py\"")
      )
    )

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedArray, Seq.empty)
    }
  }

  it should "evaluate a quote expression correctly with an array of strings" in {
    val str = """ quote(["a", "b", "c"]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedArray: WomArray = WomArray(
      Seq(
        WomString("\"a\""),
        WomString("\"b\""),
        WomString("\"c\"")
      )
    )

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedArray, Seq.empty)
    }
  }

  it should "evaluate a quote expression correctly with an array of strings that are already in quotes" in {
    val str = """ quote(["\"a\"", "\"b", "c\""]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedArray: WomArray = WomArray(
      Seq(
        WomString("\"\"a\"\""),
        WomString("\"\"b\""),
        WomString("\"c\"\"")
      )
    )

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedArray, Seq.empty)
    }
  }

  it should "evaluate a squote expression correctly with an empty array" in {
    val str = """ squote([]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedArray: WomArray = WomArray(Seq())

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedArray, Seq.empty)
    }
  }

  it should "evaluate a squote expression correctly with an array of integers" in {
    val str = """ squote([1, 2, 3]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedArray: WomArray = WomArray(
      Seq(
        WomString("\'1\'"),
        WomString("\'2\'"),
        WomString("\'3\'")
      )
    )

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedArray, Seq.empty)
    }
  }

  it should "evaluate a squote expression correctly with an array of booleans" in {
    val str = """ squote([true, false, true]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedArray: WomArray = WomArray(
      Seq(
        WomString("\'true\'"),
        WomString("\'false\'"),
        WomString("\'true\'")
      )
    )

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedArray, Seq.empty)
    }
  }

  it should "evaluate a squote expression correctly with an array of floats" in {
    val str = """ squote([1.1, 2.2, 3]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedArray: WomArray = WomArray(
      Seq(
        WomString("\'1.1\'"),
        WomString("\'2.2\'"),
        WomString("\'3.0\'")
      )
    )

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedArray, Seq.empty)
    }
  }

  it should "evaluate a squote expression correctly with an array of files" in {
    val str = """ squote(["/home/someFile.txt", "/rootFile.exe", "./anotherFile.py"]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedArray: WomArray = WomArray(
      Seq(
        WomString("\'/home/someFile.txt\'"),
        WomString("\'/rootFile.exe\'"),
        WomString("\'./anotherFile.py\'")
      )
    )

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedArray, Seq.empty)
    }
  }

  it should "evaluate a squote expression correctly with an array of strings" in {
    val str = """ squote(["a", "b", "c"]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedArray: WomArray = WomArray(
      Seq(
        WomString("\'a\'"),
        WomString("\'b\'"),
        WomString("\'c\'")
      )
    )

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedArray, Seq.empty)
    }
  }

  it should "evaluate a squote expression correctly with an array of strings that are already in quotes" in {
    val str = """ squote(["'a'", "'b", "c'"]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val expectedArray: WomArray = WomArray(
      Seq(
        WomString("""''a''"""),
        WomString("""''b'"""),
        WomString("""'c''""")
      )
    )

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedArray, Seq.empty)
    }
  }

  it should "evaluate an unzip expression correctly" in {
    val str = """ unzip([("one", 1),("two", 2),("three", 3)]) """
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val left: WomArray =
      WomArray(WomArrayType(WomStringType), Seq(WomString("one"), WomString("two"), WomString("three")))
    val right: WomArray = WomArray(WomArrayType(WomIntegerType), Seq(WomInteger(1), WomInteger(2), WomInteger(3)))
    val expectedPair: WomPair = WomPair(left, right)

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedPair, Seq.empty)
    }
  }

  it should "evaluate an unzip on an empty collection correctly" in {
    val str = """ unzip([])"""
    val expr = fromString[ExpressionElement](str, parser.parse_e)

    val left: WomArray = WomArray(WomArrayType(WomAnyType), Seq())
    val right: WomArray = WomArray(WomArrayType(WomAnyType), Seq())
    val expectedPair: WomPair = WomPair(left, right)

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedPair, Seq.empty)
    }
  }

  it should "fail to evaluate unzip on invalid pair" in {
    val invalidPair = """ unzip([()])"""
    val invalidPairExpr = fromString[ExpressionElement](invalidPair, parser.parse_e)
    invalidPairExpr.shouldBeInvalid("Failed to parse expression (reason 1 of 1): No WDL support for 0-tuples")
  }

  it should "fail to evaluate unzip on heterogeneous pairs" in {
    val invalidPair = """ unzip([ (1, 11.0), ([1,2,3], 2.0) ])"""
    val invalidPairExpr = fromString[ExpressionElement](invalidPair, parser.parse_e)
    invalidPairExpr.map(e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None)
        .shouldBeInvalid(
          "Could not construct array of type WomMaybeEmptyArrayType(WomPairType(WomIntegerType,WomFloatType)) with this value: List(WomPair(WomInteger(1),WomFloat(11.0)), WomPair([1, 2, 3],WomFloat(2.0)))"
        )
    )
  }

  it should "evaluate a struct literal" in {
    val literal = """ Animal{type: "dog", barks: false}"""
    val expectedValue = WomObject(Map("type" -> WomString("dog"), "barks" -> WomBoolean(false)))
    val expr = fromString[ExpressionElement](literal, parser.parse_e)

    expr.shouldBeValidPF { case e =>
      e.evaluateValue(Map.empty, NoIoFunctionSet, None) shouldBeValid EvaluatedValue(expectedValue, Seq.empty)
    }
  }
}
