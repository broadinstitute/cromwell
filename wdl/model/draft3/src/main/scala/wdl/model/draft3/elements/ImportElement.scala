package wdl.model.draft3.elements

final case class ImportElement(importUrl: String,
                         alias: Option[String]) extends LanguageElement
