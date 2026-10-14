version 1.2

workflow multiline_strings {
  String name = "world"

  output {
    String simple = <<<
      This is a
      multi-line string!
    >>>
    String continued = <<<
        hello  \
            world
    >>>
    String with_quotes = <<<
      multi-line string \
      with 'single' and "double" quotes
    >>>
    String with_placeholder = <<<
      hello ~{name}
    >>>
    String blank_lines = <<<

      this is a

        multi-line string

    >>>
  }
}
