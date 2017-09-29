package cwl

import wdl.types._
import wdl.values._
import wom.expression.IoFunctionSet

// TODO: WOM: We need a CWL file as defined by http://www.commonwl.org/v1.0/CommandLineTool.html#File
case class CwlFile(path: String)

object CwlFileEvaluator {

  def pathsToArrayOfCwlFileMaps(paths: Seq[String], ioFunctionSet: IoFunctionSet, loadContents: Boolean): WdlArray = {
    cwlFilesToWdlArray(paths map CwlFile.apply, ioFunctionSet, loadContents)
  }

  def cwlFilesToWdlArray(cwlFiles: Seq[CwlFile],
                         ioFunctionSet: IoFunctionSet,
                         loadContents: Boolean): WdlArray = {

    val wdlMapType = WdlMapType(WdlStringType, WdlStringType)
    val wdlMaps = cwlFiles map { cwlFile =>
      // TODO: WOM: basename/dirname/size/checksum/etc.

      val contents: Map[WdlValue, WdlValue] =
        if (loadContents) Map(WdlString("contents") -> WdlString(load64KiB(cwlFile.path, ioFunctionSet))) else Map.empty

      val wdlKeyValues: Map[WdlValue, WdlValue] = Map(
        WdlString("location") -> WdlString(cwlFile.path)
      ) ++ contents

      WdlMap(wdlMapType, wdlKeyValues)
    }

    WdlArray(WdlArrayType(wdlMapType), wdlMaps)
  }

  private def load64KiB(path: String, ioFunctionSet: IoFunctionSet): String = {
    // This suggests the IoFunctionSet should have a length-limited read API as both CWL and WDL support this concept.
    // val content = ioFunctionSet.readFile(path)
    val content = better.files.File(path).bytes.take(64 * 1024).toArray
    new String(content)
  }
}
