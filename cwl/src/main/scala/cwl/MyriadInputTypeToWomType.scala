package cwl

import shapeless.Poly1

object MyriadInputTypeToWomType extends Poly1 {
  implicit def x = at[CwlType]
    implicit def x = at[InputRecordSchema]
    implicit def x = at[InputEnumSchema]
    implicit def x = at[InputArraySchema]
    implicit def x = at[String]
  implicit def x = at[Array[]]

  /*
  CwlType :+:
    InputRecordSchema :+:
    InputEnumSchema :+:
    InputArraySchema :+:
    String :+:
    Array[
      CwlType :+:
        InputRecordSchema :+:
        InputEnumSchema :+:
        InputArraySchema :+:
        String :+:
        CNil
      ] :+:
    CNil
    */
}

object MyriadInputInnerTypeToWomType extends Poly1 {

}
