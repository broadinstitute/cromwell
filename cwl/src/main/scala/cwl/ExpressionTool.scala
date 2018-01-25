package cwl

import cwl.CwlVersion._
import cwl.ExpressionTool.{ExpressionToolInputParameter, ExpressionToolOutputParameter}
import shapeless.Coproduct

case class ExpressionTool(
                           inputs: Array[ExpressionToolInputParameter] = Array.empty,
                           outputs: Array[ExpressionToolOutputParameter] = Array.empty,
                           `class`: String,
                           expression: StringOrExpression,
                           id: Option[String] = None,
                           requirements: Option[Array[Requirement]] = None,
                           hints: Option[Array[Hint]] = None,
                           label: Option[String] = None,
                           doc: Option[String] = None,
                           cwlVersion: Option[CwlVersion] = None) {
  def asCwl = Coproduct[Cwl](this)
}

object ExpressionTool {

  case class ExpressionToolInputParameter(id: String,
                                          label: Option[String] = None,
                                          secondaryFiles: Option[SecondaryFiles] = None,
                                          format: Option[InputParameterFormat] = None,
                                          streamable: Option[Boolean] = None,
                                          doc: Option[Doc] = None,
                                          inputBinding: Option[InputCommandLineBinding] = None,
                                          default: Option[CwlAny] = None,
                                          `type`: Option[MyriadInputType] = None) extends InputParameter

  case class ExpressionToolOutputParameter(id: String,
                                           label: Option[String] = None,
                                           secondaryFiles: Option[SecondaryFiles] = None,
                                           format: Option[OutputParameterFormat] = None,
                                           streamable: Option[Boolean] = None,
                                           doc: Option[Doc] = None,
                                           outputBinding: Option[CommandOutputBinding] = None,
                                           `type`: MyriadOutputType)
}
