package cwl

import shapeless.{:+:, CNil, Coproduct}
import cwl.CwlVersion._

case class ExpressionTool(
                           inputs: Array[InputParameter] = Array.empty,
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

case class ExpressionToolOutputParameter(
                                          id: String,
                                          label: Option[String] = None,
                                          secondaryFiles: Option[SecondaryFiles] = None,
                                          format:
                                            Option[
                                              Expression :+:
                                              String :+:
                                              Array[String] :+:
                                              CNil] = None,
                                          streamable: Option[Boolean] = None,
                                          doc:
                                            Option[
                                              String :+:
                                              Array[String] :+:
                                              CNil] = None,
                                          outputBinding: Option[CommandOutputBinding] = None,
                                          `type`: MyriadOutputType)
