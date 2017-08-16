package cromwell.backend.wdl

import cromwell.core.path.Path
import wdl4s.wdl.TsvSerializable
import wdl4s.wdl.expression.WdlStandardLibraryFunctions
import wdl4s.wdl.types._
import wdl4s.wdl.values._

import scala.util.{Failure, Try}

trait WriteFunctions { this: WdlStandardLibraryFunctions =>

  /**
    * Directory that will be used to write files.
    */
  def writeDirectory: Path

  private lazy val _writeDirectory = writeDirectory.createPermissionedDirectories()

  def writeTempFile(path: String,prefix: String,suffix: String,content: String): String = throw new NotImplementedError("This method is not used anywhere and should be removed")

  private def writeContent(baseName: String, content: String): Try[WdlFile] = {
    val tmpFile = _writeDirectory / s"${baseName}_${content.md5Sum}.tmp"

    Try {
      if (tmpFile.notExists) tmpFile.write(content)
    } map { _ =>
      WdlFile(tmpFile.pathAsString)
    }
  }

  private def writeToTsv[A <: WdlValue with TsvSerializable](functionName: String, params: Seq[Try[WdlValue]], defaultIfOptionalEmpty: A): Try[WdlFile] = {
    val wdlClass = defaultIfOptionalEmpty.getClass
    def castOrDefault(wdlValue: WdlValue): A = wdlValue match {
      case WdlOptionalValue(_, None) => defaultIfOptionalEmpty
      case WdlOptionalValue(_, Some(v)) => wdlClass.cast(v)
      case _ => wdlClass.cast(wdlValue)
    }

    for {
      singleArgument <- extractSingleArgument(functionName, params)
      downcast <- Try(castOrDefault(singleArgument))
      tsvSerialized <- downcast.tsvSerialize
      file <- writeContent(functionName, tsvSerialized)
    } yield file
  }

  override def write_lines(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv("write_lines", params, WdlArray(WdlArrayType(WdlStringType), List.empty[WdlValue]))
  override def write_map(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv("write_map", params, WdlMap(WdlMapType(WdlStringType, WdlStringType), Map.empty[WdlValue, WdlValue]))
  override def write_object(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv("write_object", params, WdlObject(Map.empty[String, WdlValue]))
  override def write_objects(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv("write_objects", params, WdlArray(WdlArrayType(WdlObjectType), List.empty[WdlObject]))
  override def write_tsv(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv("write_tsv", params, WdlArray(WdlArrayType(WdlStringType), List.empty[WdlValue]))
  override def write_json(params: Seq[Try[WdlValue]]): Try[WdlFile] = Failure(new NotImplementedError(s"write_json() not implemented yet"))
}
