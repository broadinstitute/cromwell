package cromwell.engine.backend.jes

import java.nio.file.Path
import java.security.MessageDigest

import cromwell.binding.values._
import cromwell.engine.EngineFunctions
import cromwell.util.GcsUtil
import cromwell.util.google.GoogleCloudStoragePath

import org.apache.commons.codec.binary.Base64

import scala.util.{Success, Try}

/**
 * Implementation of engine functions for the JES backend.
 */
case class JesEngineFunctions(callDir: GoogleCloudStoragePath, jesConnection: JesInterface) extends EngineFunctions {

  private def readFromPath(value: String): String = {
    val tryParse = GoogleCloudStoragePath.parse(value)
    val gcsPathToUse: GoogleCloudStoragePath = if (tryParse.isSuccess) { tryParse.get } else { GoogleCloudStoragePath(gcsPathFromAnyString(value)) }
    jesConnection.storage.slurpFile(gcsPathToUse)
  }
  
  /**
   * Read the entire contents of a file from the specified `WdlValue`, where the file can be
   * specified either as a path via a `WdlString` (with magical handling of "stdout"), or
   * directly as a `WdlFile`.
   *
   * @throws UnsupportedOperationException for an unrecognized file reference, as this is intended
   *                                       to be wrapped in a `Try`.
   */
  private def fileContentsToString(value: WdlValue): String = {
    value match {
      case f: WdlFile => {
        readFromPath(f.value)
      }
      case f: WdlString => {
        readFromPath(f.value)
      }
      case e => throw new UnsupportedOperationException("Unsupported argument " + e + " (expected JES URI)")
    }
  }

  protected def read_lines(params: Seq[Try[WdlValue]]): Try[WdlArray] = ???

  /**
   * Try to read a string from the file referenced by the specified `WdlValue`.
   */
  override protected def read_string(params: Seq[Try[WdlValue]]): Try[WdlString] = {
    for {
      singleArgument <- extractSingleArgument(params)
      string = fileContentsToString(singleArgument)
    } yield WdlString(string)
  }

  /**
   * Try to read an integer from the file referenced by the specified `WdlValue`.
   */
  override protected def read_int(params: Seq[Try[WdlValue]]): Try[WdlInteger] =
    read_string(params).map { s => WdlInteger(s.value.trim.toInt) }

  override protected def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = {
    val newPath = GoogleCloudStoragePath(callDir.bucket, callDir.objectName + "/" + JesBackend.LocalStdoutValue)
    Success(WdlFile(newPath.toString))
  }

  override protected def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = {
    val newPath = GoogleCloudStoragePath(callDir.bucket, callDir.objectName + "/" + JesBackend.LocalStderrValue)
    Success(WdlFile(newPath.toString))
  }

  def gcsPathFromAnyString(value: String) = {
    callDir + md5(value)
  }

  def md5(value: String): String = {
    new String(Base64.encodeBase64(MessageDigest.getInstance("MD5").digest(value.getBytes)))
  }
}
