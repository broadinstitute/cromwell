package cwl

import cwl.CommandLineTool.CommandInputParameter
import org.scalatest.{FlatSpec, Matchers}

class CommandLineToolSpec extends FlatSpec with Matchers {

  behavior of "CommandLineTool"

  it should "filter input arguments when not bound to command line" in {

    val a = CommandInputParameter("a", inputBinding = None)
    val b = CommandInputParameter("b", inputBinding = Some(CommandLineBinding()))

    CommandLineTool.orderedForCommandLine(Array(a,b)) shouldBe Seq(b)
  }

  it should "order input arguments that are bound to command line" in {
    val a = CommandInputParameter("a", inputBinding = None)
    val b = CommandInputParameter("b", inputBinding = Some(CommandLineBinding(position = Some(1))))
    val c = CommandInputParameter("c", inputBinding = Some(CommandLineBinding(position = Some(2))))
    val d = CommandInputParameter("d", inputBinding = Some(CommandLineBinding(position = Some(3))))

    CommandLineTool.orderedForCommandLine(Array(d,a,c,b)) shouldBe Seq(b,c,d)
  }

  it should "default input arguments with Command Line Binding but no position should be assigned position 0" in {
    val a = CommandInputParameter("a", inputBinding = None)
    val b = CommandInputParameter("b", inputBinding = Some(CommandLineBinding(position = Some(1))))
    val c = CommandInputParameter("c", inputBinding = Some(CommandLineBinding()))

    CommandLineTool.orderedForCommandLine(Array(a,c,b)) shouldBe Seq(c,b)
  }
}
