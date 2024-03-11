package wdl.model.draft3.elements
import wom.SourceFileLocation

final case class CallElement(callableReference: String,
                             alias: Option[String],
                             afters: Vector[String],
                             body: Option[CallBodyElement],
                             override val sourceLocation: Option[SourceFileLocation]
) extends LanguageElement
    with WorkflowGraphElement {
  override def toString: String =
    s"""Call "$callableReference${alias.map(alias => s" as $alias").getOrElse("")}${afters
        .map(after => s" after $after")
        .mkString}""""
}
