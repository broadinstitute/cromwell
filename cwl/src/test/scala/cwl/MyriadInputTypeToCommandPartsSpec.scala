package cwl

import cats.data.NonEmptyList
import cwl.CommandLineTool.{CommandBindingSortingKey, CommandInputParameter, StringOrInt}
import cwl.SchemaDefRequirement.SchemaDefTypes
import cwl.internal.CommandPartSortingAlgorithm
import org.scalatest.{FlatSpec, Matchers}
import shapeless.Coproduct
import shapeless.syntax.inject._
import wom.types.{WomCompositeType, WomCoproductType, WomEnumerationType, WomIntegerType}
import wom.values.{WomCoproductValue, WomEnumerationValue, WomInteger, WomObject, WomValue}


class MyriadInputTypeToCommandPartsSpec extends FlatSpec with Matchers{

  "two types w/ overlapping field names" should "only produce command parts from one type" in {

    val miit = CwlType.Int.inject[MyriadInputInnerType]
    val miit2 = CwlType.Null.inject[MyriadInputInnerType]

    val mit = Array(miit, miit2).inject[MyriadInputType]

    val cip = CommandInputParameter(id = "x", `type` = Some(mit))

    val enum1 = "map1"
    val enum2 = "map2"

    val wetMap1 = WomEnumerationType(NonEmptyList.one(enum1))
    val wetMap2 = WomEnumerationType(NonEmptyList.one(enum2))


    val commonEnumFieldName = "m"
    val commonIntField = "i"

    val otherField = "o"

    val ies1 = InputEnumSchema(commonEnumFieldName, Array(enum1))

    val wit = WomIntegerType

    val wct = WomCompositeType(Map(commonEnumFieldName -> wetMap1, commonIntField -> wit, otherField -> wit))
    val wct2 = WomCompositeType(Map(commonEnumFieldName -> wetMap2, commonIntField -> wit))

    val coproductType = WomCoproductType(NonEmptyList.of(wct, wct2))

    val womIntValue = 4

    val map1 = InputRecordSchema(
      name = "s#Map1",
      fields = Some(Array(
        InputRecordField(s"x#Map1/$commonEnumFieldName", wetMap1.inject[MyriadInputInnerType].inject[MyriadInputType], inputBinding = Some(InputCommandLineBinding(position = Some(0)))),
        InputRecordField(s"x#Map1/$commonIntField", wit.inject[MyriadInputInnerType].inject[MyriadInputType], inputBinding = Some(InputCommandLineBinding(position = Some(2)))),
        InputRecordField(s"x#Map1/$otherField", wit.inject[MyriadInputInnerType].inject[MyriadInputType], inputBinding = Some(InputCommandLineBinding(position = Some(2))))
      ))).inject[SchemaDefTypes]

    val map2 = InputRecordSchema(
      name = "s#Map2",
      fields = Some(Array(
        InputRecordField(s"x#Map2/$commonEnumFieldName", wetMap2.inject[MyriadInputInnerType].inject[MyriadInputType], inputBinding = Some(InputCommandLineBinding(position = Some(0)))),
        InputRecordField(s"x#Map2/$commonIntField", wit.inject[MyriadInputInnerType].inject[MyriadInputType], inputBinding = Some(InputCommandLineBinding(position = Some(2)))),
      ))).inject[SchemaDefTypes]

    val stage = InputRecordSchema(
      name = "s#Stage",
      fields = Some(Array(
        InputRecordField(
          s"x#Stage/$commonEnumFieldName",
          Array(
            wetMap1.inject[MyriadInputInnerType],
            wetMap2.inject[MyriadInputInnerType]
          ).inject[MyriadInputType], inputBinding = Some(InputCommandLineBinding(position = Some(2))))
      ))).inject[SchemaDefTypes]

    val x: MyriadInputType =  Array(wct.inject[MyriadInputInnerType], wct2.inject[MyriadInputInnerType]).inject[MyriadInputType]

    //(Option[InputCommandLineBinding], WomValue, CommandBindingSortingKey, Boolean, ExpressionLib, SchemaDefRequirement) => CommandPartsList
    val wo = WomObject(Map(
      commonEnumFieldName -> WomEnumerationValue(wetMap1,enum1),
      commonIntField -> WomInteger(womIntValue),
      otherField -> WomInteger(womIntValue)
    ))

    val iclb: Option[InputCommandLineBinding] = Some(InputCommandLineBinding(position = Some(2)))
    val womValue: WomValue = WomCoproductValue(coproductType, wo)
    val cbsk: CommandBindingSortingKey = CommandBindingSortingKey(List(1.inject[StringOrInt]))
    val bool: Boolean = false
    val expressionLib: ExpressionLib = Vector.empty
    val sdr: SchemaDefRequirement = SchemaDefRequirement(types = Array(
      map1.inject[SchemaDefTypes],
      map2.inject[SchemaDefTypes],
      stage.inject[SchemaDefTypes]
    ))
    val result: Seq[CommandLineTool.SortKeyAndCommandPart] = x.fold(MyriadInputTypeToSortedCommandParts).apply(iclb, womValue, cbsk, bool, expressionLib, sdr)

    result.size shouldBe 3
  }

}
