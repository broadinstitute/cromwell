package cromwell.database.migration

import java.io.{ByteArrayInputStream, IOException}
import java.nio.charset.Charset
import java.util.zip.GZIPInputStream

import org.apache.commons.codec.binary.Base64
import org.apache.commons.io.IOUtils
import wdl.draft2.model.types.WdlFlavoredWomType._
import wom.types.{WomPrimitiveType, WomType}

import scala.util.Try

private[migration] object WdlTransformation {

  def inflate(value: String): Try[String] = Try {
    Option(value) match {
      case Some(v) =>
        IOUtils.toString(new GZIPInputStream(new ByteArrayInputStream(Base64.decodeBase64(v))), Charset.defaultCharset)
      case None => null
    }
  } recover { case _: IOException =>
    value
  }

  def coerceStringToWdl(wdlString: String, womType: WomType) = womType match {
    case _ if wdlString == null => null
    case p: WomPrimitiveType => p.coerceRawValue(wdlString).get
    case o => o.fromWorkflowSource(wdlString)
  }
}
