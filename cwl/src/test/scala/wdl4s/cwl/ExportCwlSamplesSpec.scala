package wdl4s.cwl

import io.circe.Printer
import org.scalatest.{FlatSpec, Matchers}
import shapeless.Coproduct

class ExportCwlSamplesSpec extends FlatSpec with Matchers {

  it should "export 1st tool" in {
    val tool =
      CommandLineTool(
        inputs = Coproduct[CommandLineTool.Inputs](Map("message" -> CommandInputParameter(
          id = None,
          label = None,
          secondaryFiles = None,
          format = None,
          streamable = None,
          doc = None,
          inputBinding = Option(CommandLineBinding(
            loadContents = None,
            position = Option(1),
            prefix = None,
            separate = None,
            itemSeparator = None,
            valueFrom = None,
            shellQuote = None)),
          default = None,
          `type` = None
        ))),
        outputs = Coproduct[CommandLineTool.Outputs](Array.empty[CommandOutputParameter]),
        `class` = "CommandLineTool",
        id = None,
        requirements = None,
        hints = None,
        label = None,
        doc = None,
        cwlVersion = Option(CwlVersion.Version1),
        baseCommand = None,
        arguments = None,
        stdin = None,
        stderr = None,
        stdout = None,
        successCodes = None,
        temporaryFailCodes = None,
        permanentFailCodes = None)
    val toolJson = encodeCwlCommandLineTool(tool)
    val printer = new Printer(preserveOrder = true, dropNullKeys = true, indent = "  ")
    val toolJsonString = printer.pretty(toolJson)
    val expectedToolJsonString =
      """{"inputs":{"message":{"inputBinding":{"position":1}}},"outputs":[],"class":"CommandLineTool","cwlVersion":"v1.0"}"""
    toolJsonString shouldBe expectedToolJsonString
  }

}
