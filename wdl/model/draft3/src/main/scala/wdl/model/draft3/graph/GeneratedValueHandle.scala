package wdl.model.draft3.graph

import wom.types.WomType

sealed trait GeneratedValueHandle {
  def linkableName: String
  def womType: WomType
}

final case class GeneratedIdentifierValueHandle(linkableName: String, womType: WomType) extends GeneratedValueHandle
final case class GeneratedCallOutputValueHandle(taskName: String, outputName: String, womType: WomType) extends GeneratedValueHandle {
  override def linkableName: String = s"$taskName.$outputName"
}
