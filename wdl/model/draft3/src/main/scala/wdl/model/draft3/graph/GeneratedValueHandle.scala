package wdl.model.draft3.graph

import wom.types.WomType

sealed trait GeneratedValueHandle {
  def linkableName: String
  def womType: WomType
}

final case class GeneratedIdentifierValueHandle(linkableName: String, womType: WomType) extends GeneratedValueHandle
final case class GeneratedCallOutputValueHandle(callName: String, outputName: String, womType: WomType) extends GeneratedValueHandle {
  override def linkableName: String = s"$callName.$outputName"
}

// A handle representing a call's completion
final case class GeneratedCallOutputAsStructHandle(finishedCallName: String, structType: WomType) extends GeneratedValueHandle {
  override def linkableName: String = s"$finishedCallName"
  override val womType: WomType = structType
}
