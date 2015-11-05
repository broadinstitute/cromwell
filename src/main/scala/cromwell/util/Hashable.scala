package cromwell.util

trait Hashable extends Any {
  def md5Sum: String
}
