package cromwell.docker

import scala.util.{Failure, Success, Try}

object DockerHashResult {
  // See https://docs.docker.com/registry/spec/api/#/content-digests
  val DigestRegex = """([a-zA-Z0-9_+.-]+):([a-zA-Z0-9]+)""".r
  
  def fromString(str: String): Try[DockerHashResult] = {
    str match {
      case DigestRegex(alg, hash) => Success(DockerHashResult(alg, hash))
      case _ => Failure(new IllegalArgumentException(s"Hash value $str does not have the expected 'algorithm:hash' syntax"))
    }
  }
}

case class DockerHashResult(hashAlgorithm: String, hashValue: String) {
  val algorithmAndHash = s"$hashAlgorithm:$hashValue"
}
