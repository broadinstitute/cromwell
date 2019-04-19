package wom.values

import org.scalatest.{FlatSpec, Matchers}
import wom.types.{WomIntegerType, WomOptionalType}

class WomOptionalValueSpec extends FlatSpec with Matchers {

  behavior of "WomOptionalValue flattening"

  // Unnested optionals flatten to themseles:

  it should "flatten an unnested Some in an Int? to itself" in {
    val opt = WomOptionalValue(WomInteger(75))
    opt.womType.stableName should be("Int?")

    val flattened = opt.flattenOptional
    flattened should be(WomOptionalValue(WomInteger(75)))
    flattened.womType.stableName should be("Int?")
  }

  it should "flatten an unnested None in an Int? to itself" in {
    val opt = WomOptionalValue(WomIntegerType, None)
    opt.womType.stableName should be("Int?")

    val flattened = opt.flattenOptional
    flattened should be(WomOptionalValue(WomIntegerType, None))
    flattened.womType.stableName should be("Int?")
  }

  // Once nested optionals:

  it should "flatten an once-nested Some in an Int?? to an unnested Some in an Int?" in {
    val opt = WomOptionalValue(WomOptionalValue(WomInteger(75)))
    opt.womType.stableName should be("Int??")

    val flattened = opt.flattenOptional
    flattened should be(WomOptionalValue(WomInteger(75)))
    flattened.womType.stableName should be("Int?")
  }

  it should "flatten an unnested None in an Int?? to an unnested None in an Int?" in {
    val opt = WomOptionalValue(WomOptionalType(WomIntegerType), None)
    opt.womType.stableName should be("Int??")

    val flattened = opt.flattenOptional
    flattened should be(WomOptionalValue(WomIntegerType, None))
    flattened.womType.stableName should be("Int?")
  }

  it should "flatten an once-nested None in an Int?? to an unnested None in an Int?" in {
    val opt = WomOptionalValue(WomOptionalValue(WomIntegerType, None))
    opt.womType.stableName should be("Int??")

    val flattened = opt.flattenOptional
    flattened should be(WomOptionalValue(WomIntegerType, None))
    flattened.womType.stableName should be("Int?")
  }

  // Twice-nested optionals:

  it should "flatten a twice-nested Some in an Int??? to an unnested Some in an Int?" in {
    val opt = WomOptionalValue(WomOptionalValue(WomOptionalValue(WomInteger(75))))
    opt.womType.stableName should be("Int???")

    val flattened = opt.flattenOptional
    flattened should be(WomOptionalValue(WomInteger(75)))
    flattened.womType.stableName should be("Int?")
  }

  it should "flatten an unnested None in an Int??? to an unnested None in an Int?" in {
    val opt = WomOptionalValue(WomOptionalType(WomOptionalType(WomIntegerType)), None)
    opt.womType.stableName should be("Int???")

    val flattened = opt.flattenOptional
    flattened should be(WomOptionalValue(WomIntegerType, None))
    flattened.womType.stableName should be("Int?")
  }

  it should "flatten an once-nested None in an Int??? to an unnested None in an Int?" in {
    val opt = WomOptionalValue(WomOptionalValue(WomOptionalType(WomIntegerType), None))
    opt.womType.stableName should be("Int???")

    val flattened = opt.flattenOptional
    flattened should be(WomOptionalValue(WomIntegerType, None))
    flattened.womType.stableName should be("Int?")
  }

  it should "flatten a twice-nested None in an Int??? to an unnested None in an Int?" in {
    val opt = WomOptionalValue(WomOptionalValue(WomOptionalValue(WomIntegerType, None)))
    opt.womType.stableName should be("Int???")

    val flattened = opt.flattenOptional
    flattened should be(WomOptionalValue(WomIntegerType, None))
    flattened.womType.stableName should be("Int?")
  }
}
