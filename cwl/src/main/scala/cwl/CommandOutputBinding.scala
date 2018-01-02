package cwl

import cwl.CommandOutputBinding.Glob
import shapeless.{:+:, CNil}

/** @see <a href="http://www.commonwl.org/v1.0/Workflow.html#CommandOutputBinding">CommandOutputBinding</a> */
case class CommandOutputBinding(
                                 glob: Option[Glob] = None,
                                 loadContents: Option[Boolean] = None,
                                 outputEval: Option[StringOrExpression] = None)

object CommandOutputBinding {
  type Glob = Expression :+: String :+: Array[String] :+: CNil

}
