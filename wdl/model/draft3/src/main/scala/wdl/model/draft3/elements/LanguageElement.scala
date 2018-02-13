package wdl.model.draft3.elements

trait LanguageElement {
  def children: Seq[LanguageElement]
}
