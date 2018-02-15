package wdl.model.draft3.elements

case class ImportElement(importUrl: String,
                         alias: Option[String]) extends LanguageElement
