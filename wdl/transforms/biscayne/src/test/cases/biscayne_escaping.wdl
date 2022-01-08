version development

workflow escapes {
    String backslash = " \\ "
    String n = " \n "
    String t = " \t "

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

    String octal_hello = "\150\145\154\154\157"
    String hex_hello = "\x68\x65\x6C\x6c\x6F"
    String unicode_hello = "\u0068\U00000065\u006C\U0000006C\u006F"
}
