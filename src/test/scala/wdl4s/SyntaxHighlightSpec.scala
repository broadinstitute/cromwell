package wdl4s

import wdl4s.formatter.{AnsiSyntaxHighlighter, HtmlSyntaxHighlighter, SyntaxFormatter}
import org.scalatest.{Matchers, WordSpecLike}

class SyntaxHighlightSpec extends Matchers with WordSpecLike {

  "SyntaxFormatter for simple workflow" should {
    val namespace = WdlNamespace.load(
      """task t {
        |  String f
        |  Int p
        |  command {
        |    ./cmd ${f} ${p}
        |  }
        |}
        |workflow w {
        |  Array[String] a = ["foo", "bar", "baz"]
        |  call t
        |  call t as u {
        |    input: f="abc", p=2
        |  }
        |  scatter (b in a) {
        |    call t as v {
        |      input: f=t, p=3
        |    }
        |  }
        |}
      """.stripMargin
    )

    val console =
      """\u001b[38;5;214mtask\u001b[0m \u001b[38;5;253mt\u001b[0m {
        |  \u001b[38;5;33mString\u001b[0m \u001b[38;5;112mf\u001b[0m
        |  \u001b[38;5;33mInt\u001b[0m \u001b[38;5;112mp\u001b[0m
        |  \u001b[38;5;214mcommand\u001b[0m {
        |    ./cmd ${f} ${p}
        |  }
        |}
        |
        |\u001b[38;5;214mworkflow\u001b[0m \u001b[38;5;253mw\u001b[0m {
        |  \u001b[38;5;33mArray[String]\u001b[0m \u001b[38;5;112ma\u001b[0m = ["foo","bar","baz"]
        |  \u001b[38;5;214mcall\u001b[0m \u001b[38;5;253mt\u001b[0m
        |  \u001b[38;5;214mcall\u001b[0m \u001b[38;5;253mt\u001b[0m as u {
        |    input: f="abc", p=2
        |  }
        |  \u001b[38;5;214mscatter\u001b[0m (b in a) {
        |    \u001b[38;5;214mcall\u001b[0m \u001b[38;5;253mt\u001b[0m as v {
        |      input: f=t, p=3
        |    }
        |  }
        |}""".stripMargin

    val html =
      """<span class="keyword">task</span> <span class="name">t</span> {
        |  <span class="type">String</span> <span class="variable">f</span>
        |  <span class="type">Int</span> <span class="variable">p</span>
        |  <span class="section">command</span> {
        |    <span class="command">./cmd ${f} ${p}</span>
        |  }
        |}
        |
        |<span class="keyword">workflow</span> <span class="name">w</span> {
        |  <span class="type">Array[String]</span> <span class="variable">a</span> = ["foo","bar","baz"]
        |  <span class="keyword">call</span> <span class="name">t</span>
        |  <span class="keyword">call</span> <span class="name">t</span> as <span class="alias">u</span> {
        |    input: f="abc", p=2
        |  }
        |  <span class="keyword">scatter</span> (b in a) {
        |    <span class="keyword">call</span> <span class="name">t</span> as <span class="alias">v</span> {
        |      input: f=t, p=3
        |    }
        |  }
        |}""".stripMargin

    "format to console properly" in {
      val actual = new SyntaxFormatter(AnsiSyntaxHighlighter).format(namespace)
      actual shouldEqual console
    }

    "format to HTML properly" in {
      val actual = new SyntaxFormatter(HtmlSyntaxHighlighter).format(namespace)
      actual shouldEqual html
    }
  }

  "SyntaxFormatter for more feature-rich workflow" should {
    val namespace = WdlNamespace.load(
      """import "foo.wdl" as foo_ns
        |
        |task t {
        |  String f
        |  Int p
        |  command {
        |    ./cmd ${f} ${p}
        |  }
        |}
        |
        |task s {
        |  Array[File] input_file
        |  command <<<
        |    cat ${sep=' ' input_file} | awk '{s+=$1} END {print s}'
        |  >>>
        |  output {
        |    String s = read_string(stdout())
        |  }
        |}
        |
        |task r {
        |  command { python -c "import random; print(random.randint(1,100))" }
        |}
        |
        |workflow w {
        |  Int p = 2+2
        |  call t
        |  call t as u {
        |    input: f="abc", p=p
        |  }
        |}""".stripMargin,
      importResolver = (s:String) => ""
    )

    val console =
      """\u001b[38;5;214mimport\u001b[0m 'foo.wdl' as foo_ns
        |
        |\u001b[38;5;214mtask\u001b[0m \u001b[38;5;253mt\u001b[0m {
        |  \u001b[38;5;33mString\u001b[0m \u001b[38;5;112mf\u001b[0m
        |  \u001b[38;5;33mInt\u001b[0m \u001b[38;5;112mp\u001b[0m
        |  \u001b[38;5;214mcommand\u001b[0m {
        |    ./cmd ${f} ${p}
        |  }
        |}
        |
        |\u001b[38;5;214mtask\u001b[0m \u001b[38;5;253ms\u001b[0m {
        |  \u001b[38;5;33mArray[File]\u001b[0m \u001b[38;5;112minput_file\u001b[0m
        |  \u001b[38;5;214mcommand\u001b[0m <<<
        |    cat ${sep=" " input_file} | awk '{s+=$1} END {print s}'
        |  >>>
        |  \u001b[38;5;214moutput\u001b[0m {
        |    \u001b[38;5;33mString\u001b[0m \u001b[38;5;112ms\u001b[0m = \u001b[38;5;13mread_string\u001b[0m(\u001b[38;5;13mstdout\u001b[0m())
        |  }
        |}
        |
        |\u001b[38;5;214mtask\u001b[0m \u001b[38;5;253mr\u001b[0m {
        |  \u001b[38;5;214mcommand\u001b[0m {
        |    python -c "import random; print(random.randint(1,100))"
        |  }
        |}
        |
        |\u001b[38;5;214mworkflow\u001b[0m \u001b[38;5;253mw\u001b[0m {
        |  \u001b[38;5;33mInt\u001b[0m \u001b[38;5;112mp\u001b[0m = 2 + 2
        |  \u001b[38;5;214mcall\u001b[0m \u001b[38;5;253mt\u001b[0m
        |  \u001b[38;5;214mcall\u001b[0m \u001b[38;5;253mt\u001b[0m as u {
        |    input: f="abc", p=p
        |  }
        |}""".stripMargin

    val html =
      """<span class="keyword">import</span> 'foo.wdl' as foo_ns
        |
        |<span class="keyword">task</span> <span class="name">t</span> {
        |  <span class="type">String</span> <span class="variable">f</span>
        |  <span class="type">Int</span> <span class="variable">p</span>
        |  <span class="section">command</span> {
        |    <span class="command">./cmd ${f} ${p}</span>
        |  }
        |}
        |
        |<span class="keyword">task</span> <span class="name">s</span> {
        |  <span class="type">Array[File]</span> <span class="variable">input_file</span>
        |  <span class="section">command</span> <<<
        |    <span class="command">cat ${sep=" " input_file} | awk '{s+=$1} END {print s}'</span>
        |  >>>
        |  <span class="section">output</span> {
        |    <span class="type">String</span> <span class="variable">s</span> = <span class="function">read_string</span>(<span class="function">stdout</span>())
        |  }
        |}
        |
        |<span class="keyword">task</span> <span class="name">r</span> {
        |  <span class="section">command</span> {
        |    <span class="command">python -c "import random; print(random.randint(1,100))"</span>
        |  }
        |}
        |
        |<span class="keyword">workflow</span> <span class="name">w</span> {
        |  <span class="type">Int</span> <span class="variable">p</span> = 2 + 2
        |  <span class="keyword">call</span> <span class="name">t</span>
        |  <span class="keyword">call</span> <span class="name">t</span> as <span class="alias">u</span> {
        |    input: f="abc", p=p
        |  }
        |}""".stripMargin

    "format to console properly" in {
      val actual = new SyntaxFormatter(AnsiSyntaxHighlighter).format(namespace)
      actual shouldEqual console
    }

    "format to HTML properly" in {
      val actual = new SyntaxFormatter(HtmlSyntaxHighlighter).format(namespace)
      actual shouldEqual html
    }
  }
}
