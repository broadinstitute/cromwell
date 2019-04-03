version 1.0

workflow escapes {
    String backslash = "\\.gz$"
    String n = "\\n"
    String r = "\\r"
    String b = "\\b"
    String t = "\\t"
    String f = "\\f"
    String a = "\\a"
    String v = "\\v"

    String q1 = "leading text \" trailing text"
    String q2 =              "\""
    String q3 =            "  \"  "
    String q4 = "leading text \' trailing text"
    String q5 =              "\'"
    String q6 =            "  \'  "

    String sq1 = 'leading text \" trailing text'
    String sq2 =              '\"'
    String sq3 =            '  \"  '
    String sq4 = 'leading text \' trailing text'
    String sq5 =              '\''
    String sq6 =            '  \'  '
}
