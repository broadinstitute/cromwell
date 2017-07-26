package wdl4s.cwl

import shapeless.{:+:, CNil}
import wdl4s.cwl.CwlVersion._

case class ExpressionTool(
                           inputs: Array[InputParameter] = Array.empty,
                           outputs: Array[ExpressionToolOutputParameter] = Array.empty,
                           `class`: String,
                           expression:
                             ECMAScriptExpression :+:
                             String :+:
                             CNil,
                           id: Option[String] = None,
                           requirements: Option[Array[Requirement]] = None,
                           hints: Option[Array[CwlAny]] = None, //TODO should be Any
                           label: Option[String] = None,
                           doc: Option[String] = None,
                           cwlVersion: Option[CwlVersion] = None
)

case class ExpressionToolOutputParameter(
                                          id: String,
                                          label: Option[String] = None,
                                          secondaryFiles:
                                            Option[ECMAScriptExpression :+:
                                            String :+:
                                            Array[
                                              ECMAScriptExpression :+:
                                              String :+:
                                              CNil] :+:
                                            CNil] = None,
                                          format:
                                            Option[
                                              ECMAScriptExpression :+:
                                              String :+:
                                              Array[String] :+:
                                              CNil] = None,
                                          streamable: Option[Boolean],
                                          doc:
                                            Option[
                                              String :+:
                                              Array[String] :+:
                                              CNil] = None,
                                          outputBinding: Option[CommandOutputBinding] = None,
                                          `type`: MyriadOutputType
) {
  type Id = String
  type `type` = MyriadOutputType
}
