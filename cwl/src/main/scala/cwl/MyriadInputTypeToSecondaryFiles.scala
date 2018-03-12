package cwl

import cwl.CwlType.CwlType
import shapeless.Poly1

// IAS.secondaryFiles are NOT listed in 1.0 spec, but according to jgentry they will be, maybe
object MyriadInputTypeToSecondaryFiles extends Poly1 {
  implicit val caseMyriadInputInnerType: Case.Aux[MyriadInputInnerType, Option[SecondaryFiles]] = at {
    _.fold(MyriadInputInnerTypeToSecondaryFiles)
  }

  implicit val caseArrayMyriadInputInnerType: Case.Aux[Array[MyriadInputInnerType], Option[SecondaryFiles]] = at {
    _.toStream.flatMap(_.fold(MyriadInputInnerTypeToSecondaryFiles)).headOption
  }
}

object MyriadInputInnerTypeToSecondaryFiles extends Poly1 {
  implicit val caseCwlType: Case.Aux[CwlType, Option[SecondaryFiles]] = at { _ => None }
  implicit val caseInputRecordSchema: Case.Aux[InputRecordSchema, Option[SecondaryFiles]] = at { _ => None }
  implicit val caseInputEnumSchema: Case.Aux[InputEnumSchema, Option[SecondaryFiles]] = at { _ => None }
  implicit val caseInputArraySchema: Case.Aux[InputArraySchema, Option[SecondaryFiles]] = at {
    inputArraySchema => inputArraySchema.secondaryFiles
  }
  implicit val caseString: Case.Aux[String, Option[SecondaryFiles]] = at { _ => None }
}
