package cromwell.engine.backend

import wom.graph.GraphNodePort.OutputPort
import wom.types._
import wom.values._

package object mock {

  // This is used by stubbed backends that are to be used in tests to prepare dummy outputs for job
  def taskOutputToJobOutput(taskOutput: OutputPort) = taskOutput -> sampleValue(taskOutput.womType)

  private def sampleValue(womType: WomType): WomValue = womType match {
    case WomIntegerType => WomInteger(3)
    case WomFloatType => WomFloat(55.55)
    case WomStringType => WomString("The rain in Spain falls mainly in the plain")
    case WomBooleanType => WomBoolean(true)
    case WomSingleFileType => WomSingleFile("/root/of/all/evil")
    case WomArrayType(memberType) => WomArray(WomArrayType(memberType), List(sampleValue(memberType)))
    case WomObjectType => WomObject(Map("a" -> WomString("1"), "b" -> WomString("2")))
    case WomMapType(keyType, valueType) =>
      WomMap(WomMapType(keyType, valueType), Map(sampleValue(keyType) -> sampleValue(valueType)))
  }
}
