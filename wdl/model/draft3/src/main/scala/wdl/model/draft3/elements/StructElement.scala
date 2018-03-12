package wdl.model.draft3.elements

final case class StructElement(name: String, entries: Seq[StructEntryElement]) extends FileBodyElement
final case class StructEntryElement(identifier: String, typeElement: TypeElement)
