package wdl.model.draft3.elements

final case class CommandSectionElement(parts: Seq[CommandSectionLine]) extends TaskSectionElement

final case class CommandSectionLine(parts: Seq[CommandPartElement])
