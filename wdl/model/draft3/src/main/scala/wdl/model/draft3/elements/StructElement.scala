package wdl.model.draft3.elements

final case class StructElement(entries: Seq[StructEntryElement]) extends FileBodyElement
final case class StructEntryElement(identifier: String, typeElement: TypeElement)
