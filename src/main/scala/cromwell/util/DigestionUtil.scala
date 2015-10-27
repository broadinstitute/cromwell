package cromwell.util

import java.security.MessageDigest
import javax.xml.bind.annotation.adapters.HexBinaryAdapter

object DigestionUtil {
  val hexBinaryAdapter = new HexBinaryAdapter()
  
  def md5Sum(input: String): String = hexBinaryAdapter.marshal(MessageDigest.getInstance("MD5").digest(input.getBytes()))
}
