package cromwell.backend.wdl

import java.nio.file.Path

import wdl4s.TsvSerializable
import wdl4s.expression.WdlStandardLibraryFunctions
import wdl4s.values._

import scala.util.{Failure, Try}

trait WriteFunctions { this: WdlStandardLibraryFunctions =>
  import better.files._

  /**
    * Directory that will be used to write files.
    */
  def writeDirectory: Path

  private lazy val _writeDirectory = File(writeDirectory).createDirectories()

  def writeTempFile(path: String,prefix: String,suffix: String,content: String): String = throw new NotImplementedError("This method is not used anywhere and should be removed")

  private def writeContent(baseName: String, content: String): Try[WdlFile] = {
    val tmpFile = _writeDirectory / s"$baseName-${content.md5Sum}.tmp"

    Try {
      if (tmpFile.notExists) tmpFile.write(content)
    } map { _ =>
      WdlFile(tmpFile.uri.toString)
    }
  }

  private def writeToTsv(functionName: String, params: Seq[Try[WdlValue]], wdlClass: Class[_ <: WdlValue with TsvSerializable]) = {
    for {
      singleArgument <- extractSingleArgument(functionName, params)
      downcast <- Try(wdlClass.cast(singleArgument))
      tsvSerialized <- downcast.tsvSerialize
      file <- writeContent(wdlClass.getSimpleName.toLowerCase, tsvSerialized)
    } yield file
  }

  override def write_lines(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv("write_lines", params, classOf[WdlArray])
  override def write_map(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv("write_map", params, classOf[WdlMap])
  override def write_object(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv("write_object", params, classOf[WdlObject])
  override def write_objects(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv("write_objects", params, classOf[WdlArray])
  override def write_tsv(params: Seq[Try[WdlValue]]): Try[WdlFile] = writeToTsv("write_tsv", params, classOf[WdlArray])
  override def write_json(params: Seq[Try[WdlValue]]): Try[WdlFile] = Failure(new NotImplementedError(s"write_json() not implemented yet"))
}
