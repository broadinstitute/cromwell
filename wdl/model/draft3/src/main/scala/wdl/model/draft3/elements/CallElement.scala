package wdl.model.draft3.elements
import wom.SourceFileLocation

final case class CallElement(callableReference: String,
                             alias: Option[String],
                             afters: Vector[String],
                             body: Option[CallBodyElement],
                             override val sourceLocation : Option[SourceFileLocation])
    extends LanguageElement with WorkflowGraphElement
