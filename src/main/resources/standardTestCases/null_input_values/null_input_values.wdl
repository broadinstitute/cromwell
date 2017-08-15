workflow null_input_values {
    String? maybeString
    Int? maybeInt
    Float? maybeFloat
    Boolean? maybeBoolean
    Array[Int?] arrayOfMaybeInts
    Array[Int]? maybeArrayOfInts
    Map[String, Int?] mapOfStringToMaybeInt
    Map[String, Int]? maybeMapOfStringToInt
    Pair[Int?, String?] pairOfMaybeIntMaybeString
    
    output {
        String? o1 = maybeString
        Int? o2 = maybeInt
        Float? o3 = maybeFloat
        Boolean? o4 = maybeBoolean
        Array[Int?] o5 = arrayOfMaybeInts
        Array[Int]? o6 = maybeArrayOfInts
        Map[String, Int?] o7 = mapOfStringToMaybeInt
        Map[String, Int]? o8 = maybeMapOfStringToInt
        Pair[Int?, String?] o9 = pairOfMaybeIntMaybeString
    }
}
