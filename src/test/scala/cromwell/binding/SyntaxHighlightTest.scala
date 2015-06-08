package cromwell.binding

import cromwell.binding.formatter.{HtmlSyntaxHighlighter, TerminalSyntaxHighlighter, SyntaxFormatter}
import org.scalatest.{FlatSpec, Matchers}

class SyntaxHighlightTest extends FlatSpec with Matchers {
  val binding = WdlBinding.process(
    s"""task t {command{./cmd $${f} $${Int p}}}
       |workflow w {
       |  call t
       |  call t as u {
       |    input: f="abc", p=2
       |  }
       |}
     """.stripMargin
  )

  "SyntaxFormatter" should "produce tagged HTML" in {
    val formatter = new SyntaxFormatter(HtmlSyntaxHighlighter)
    val actual = formatter.format(binding)
    val expected = s"""<span class="keyword">task</span> <span class="name">t</span> {
      |  <span class="section">command</span> {
      |    <span class="command">./cmd $${<span class="type">String</span> <span class="variable">f</span>} $${<span class="type">Int</span> <span class="variable">p</span>}</span>
      |  }
      |}
      |
      |<span class="keyword">workflow</span> <span class="name">w</span> {
      |  <span class="keyword">call</span> <span class="name">t</span>
      |  <span class="keyword">call</span> <span class="name">t</span> as <span class="alias">u</span> {
      |    input: f="abc", p=2
      |  }
      |}""".stripMargin
    actual shouldEqual expected
  }

  it should "produce ANSI terminal output" in {
    val formatter = new SyntaxFormatter(TerminalSyntaxHighlighter)
    val actual = formatter.format(binding)
    val expected = s"""\u001b[38;5;214mtask\u001b[0m \u001b[38;5;253mt\u001b[0m {
      |  command {
      |    ./cmd $${\u001b[38;5;33mString\u001b[0m \u001b[38;5;112mf\u001b[0m} $${\u001b[38;5;33mInt\u001b[0m \u001b[38;5;112mp\u001b[0m}
      |  }
      |}
      |
      |\u001b[38;5;214mworkflow\u001b[0m \u001b[38;5;253mw\u001b[0m {
      |  \u001b[38;5;214mcall\u001b[0m \u001b[38;5;253mt\u001b[0m
      |  \u001b[38;5;214mcall\u001b[0m \u001b[38;5;253mt\u001b[0m as u {
      |    input: f="abc", p=2
      |  }
      |}""".stripMargin
    actual shouldEqual expected
  }
}
