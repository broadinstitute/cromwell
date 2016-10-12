package cromwell.backend.impl.htcondor.caching.localization

import java.nio.file.{Files, Path, Paths}

import better.files.File
import cromwell.core.{JobOutput, _}
import wdl4s.types.{WdlArrayType, WdlFileType}
import wdl4s.values.{WdlArray, WdlSingleFile, WdlValue}

trait CachedResultLocalization {
  private[localization] def localizePathViaSymbolicLink(originalPath: Path, executionPath: Path): Path = {
    if (File(originalPath).isDirectory) throw new UnsupportedOperationException("Cannot localize directory with symbolic links.")
    else {
      File(executionPath).parent.createDirectories()
      Files.createSymbolicLink(executionPath, originalPath.toAbsolutePath)
    }
  }

  private[localization] def localizeCachedFile(executionPath: Path, output: WdlValue): WdlSingleFile = {
    val origPath = Paths.get(output.valueString)
    val newPath = executionPath.toAbsolutePath.resolve(origPath.getFileName)
    val slPath = localizePathViaSymbolicLink(origPath, newPath)
    WdlSingleFile(slPath.toString)
  }

  def localizeCachedOutputs(executionPath: Path, outputs: CallOutputs): CallOutputs = {
    outputs map { case (lqn, jobOutput) =>
      jobOutput.wdlValue.wdlType match {
        case WdlFileType => (lqn -> JobOutput(localizeCachedFile(executionPath, jobOutput.wdlValue)))
        case WdlArrayType(WdlFileType) =>
          val newArray: Seq[WdlSingleFile] = jobOutput.wdlValue.asInstanceOf[WdlArray].value map {
            localizeCachedFile(executionPath, _)
          }
          (lqn -> JobOutput(WdlArray(WdlArrayType(WdlFileType), newArray)))
        case _ => (lqn, jobOutput)
      }
    }
  }
}
