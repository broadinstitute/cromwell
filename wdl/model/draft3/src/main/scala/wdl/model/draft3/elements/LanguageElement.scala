package wdl.model.draft3.elements
import wom.SourceFileLocation

trait LanguageElement {
    val srcLoc : Option[SourceFileLocation] = None
}
