version 1.0

# This WDL contains expressions that won't be linkable, but it should be parseable and wdlom-able:
workflow input_expressions {
  input {
    Int ten = 5 + 5
    Int zero = 5 - 5
    Int twentyfive = 5 * 5
    Int one = 5 / 5
    Int zeroagain = 5 % 5

    Int tenagain = ten
    Int pair_expression_member_access = (5, "hello").left

    # Also checks that the one-item tuple is simplified out during the conversion to WDLOM:
    Boolean unary_expressions = !(+5 < -5)

    Boolean comparisons = a == b || c != d && e < f && e <= f || g > h && g >= h

    Int variableLookupMemberAccesses = a.b.c.d.e.f.g
    Int expressionMemberAccesses = (5, "hello").a.b.c.d.e.f.g

    Array[Int] is = [5, 5 * 5, -5]

    Object object_literal = object {
      a: 5,
      b: "hello",
      c: (5, "hello")
    }

    Map[String, Int] map_literal = {
      "a": 5,
      "b": 5 + 5,
      "c": 5 - 5
    }

    Int ternaryIf = if true then 5 else 5

    # 0 and 1 param functions:
    String string_read = read_string(stdout())

    # 2 param function (it's not really an Int but we're not going to validate that for this test):
    Int zipped = zip(is, is)

    # 3 param function
    String subbed = sub("hello", "hello", "hello")
  }
}
