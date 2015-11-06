package cromwell.util

import java.security.MessageDigest
import javax.xml.bind.annotation.adapters.HexBinaryAdapter

import cromwell.engine.Hash

object StringDigestion {
  private val hexBinaryAdapter = new HexBinaryAdapter()

  implicit class StringDigester(val string: String) extends AnyVal {
    def md5Sum: String = hexBinaryAdapter.marshal(MessageDigest.getInstance("MD5").digest(string.getBytes))
  }
}