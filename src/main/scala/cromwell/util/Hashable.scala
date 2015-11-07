package cromwell.util

import java.io.{File, FileInputStream}

trait Hashable extends Any {
  def md5Sum: String
}
