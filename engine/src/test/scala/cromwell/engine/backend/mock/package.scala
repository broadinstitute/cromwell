package cromwell.engine.backend

import wom.JobOutput
import wom.callable.Callable.OutputDefinition
import wom.types._
import wom.values._

package object mock {

  // This is used by stubbed backends that are to be used in tests to prepare dummy outputs for job
  def taskOutputToJobOutput(taskOutput: OutputDefinition) =
    taskOutput.name -> JobOutput(sampleValue(taskOutput.womType))

  private def sampleValue(wdlType: WdlType): WdlValue = wdlType match {
    case WdlIntegerType => WdlInteger(3)
    case WdlFloatType => WdlFloat(55.55)
    case WdlStringType => WdlString("The rain in Spain falls mainly in the plain")
    case WdlBooleanType => WdlBoolean(true)
    case WdlFileType => WdlFile("/root/of/all/evil")
    case WdlArrayType(memberType) => WdlArray(WdlArrayType(memberType), List(sampleValue(memberType)))
    case WdlObjectType => WdlObject(Map("a" -> WdlString("1"), "b" -> WdlString("2")))
    case WdlMapType(keyType, valueType) => WdlMap(WdlMapType(keyType, valueType), Map(sampleValue(keyType) -> sampleValue(valueType)))
  }
}
