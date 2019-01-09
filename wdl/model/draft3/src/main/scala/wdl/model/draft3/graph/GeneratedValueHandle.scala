package wdl.model.draft3.graph

import wom.types.{WomNothingType, WomType}

sealed trait GeneratedValueHandle {
  def linkableName: String
  def womType: WomType
}

final case class GeneratedIdentifierValueHandle(linkableName: String, womType: WomType) extends GeneratedValueHandle
final case class GeneratedCallOutputValueHandle(callName: String, outputName: String, womType: WomType) extends GeneratedValueHandle {
  override def linkableName: String = s"$callName.$outputName"
}

// A handle representing a call's completion
final case class GeneratedCallFinishedHandle(finishedCallName: String) extends GeneratedValueHandle {
  override def linkableName: String = s"$finishedCallName.__after"
  override val womType: WomType = WomNothingType
}

object GeneratedValueHandle {

  private val DotRegex = "(.*)\\.(.*)".r

  def apply(name: String, womType: WomType): GeneratedValueHandle = name match {
    case DotRegex(first, second) => GeneratedCallOutputValueHandle(first, second, womType)
    case _ => GeneratedIdentifierValueHandle(name, womType)
  }
}
