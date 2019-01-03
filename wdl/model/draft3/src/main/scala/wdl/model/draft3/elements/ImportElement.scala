package wdl.model.draft3.elements

final case class ImportElement(importUrl: String,
                               namespace: Option[String],
                               structRenames: Map[String, String]) extends LanguageElement
