package wdl4s.cwl

import shapeless.{:+:, CNil}
import wdl4s.cwl.CwlVersion._

case class ExpressionTool(
                           inputs: Array[InputParameter] :+: Map[InputParameter#Id, InputParameter#`type`] :+: Map[InputParameter#Id, InputParameter] :+: CNil,
                           outputs: Array[ExpressionToolOutputParameter] :+: Map[ExpressionToolOutputParameter#Id, ExpressionToolOutputParameter#`type`] :+: Map[ExpressionToolOutputParameter#Id, ExpressionToolOutputParameter] :+: CNil,
                           `class`: String,
                           expression: ECMAScriptExpression :+: String :+: CNil,
                           id: Option[String],
                           requirements: Option[Array[Requirement]],
                           hints: Option[Array[String]], //TODO should be Any
                           label: Option[String],
                           doc: Option[String],
                           cwlVersion: Option[CWLVersion]
)

case class ExpressionToolOutputParameter(
                                          id: String,
                                          label: Option[String],
                                          secondaryFiles: Option[ECMAScriptExpression :+: String :+: Array[ECMAScriptExpression :+: String :+: CNil] :+: CNil], //TODO Use Either if/when it works
                                          format: Option[ECMAScriptExpression :+: String :+: Array[String] :+: CNil],
                                          streamable: Option[Boolean],
                                          doc: Option[String :+: Array[String] :+: CNil],
                                          outputBinding: Option[CommandOutputBinding],
                                          `type`: MyriadOutputType
) {
  type Id = String
  type `type` = MyriadOutputType
}
