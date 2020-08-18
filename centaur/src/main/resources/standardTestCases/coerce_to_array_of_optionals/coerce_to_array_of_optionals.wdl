version 1.0

workflow coerce_to_array_of_optionals {
    # Technically, select_first and select_all take Array[String?] and we're providing Array[String], but we should be able to
    # coerce between those types:
    String first = select_first(["hello", "world"])
    Array[String] all = select_all(["hello", "world"])

    output {
        String first_out = first
        Array[String] all_out = all
    }
}
