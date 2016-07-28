package cromwell.database.migration

import java.io.{ByteArrayInputStream, IOException}
import java.util.zip.GZIPInputStream

import org.apache.commons.codec.binary.Base64
import org.apache.commons.io.IOUtils
import wdl4s.types.{WdlPrimitiveType, WdlType}

import scala.language.postfixOps
import scala.util.Try

private [migration] object WdlTransformation {

  def inflate(value: String): Try[String] = {
    Try {
      IOUtils.toString(new GZIPInputStream(new ByteArrayInputStream(Base64.decodeBase64(value))))
    } recover {
      case e: IOException => value
    }
  }

  def coerceStringToWdl(wdlString: String, wdlType: WdlType) = wdlType match {
    case p: WdlPrimitiveType => p.coerceRawValue(wdlString).get
    case o => o.fromWdlString(wdlString)
  }
}
