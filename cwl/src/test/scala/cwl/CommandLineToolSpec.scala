package cwl

import cwl.CommandLineTool.CommandInputParameter
import org.scalatest.{FlatSpec, Matchers}

class CommandLineToolSpec extends FlatSpec with Matchers {

  behavior of "CommandLineTool"

  it should "filter input arguments when not bound to command line" in {

    val a = CommandInputParameter("a", inputBinding = None)
    val b = CommandInputParameter("b", inputBinding = Some(CommandLineBinding()))

    val ordered: Seq[CommandInputParameter] = CommandLineTool.orderedForCommandLine(Array(a,b))

    ordered.contains(a) shouldBe false

    ordered.contains(b) shouldBe true
  }
}
