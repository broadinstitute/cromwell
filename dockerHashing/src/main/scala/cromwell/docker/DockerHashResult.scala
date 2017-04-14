package cromwell.docker

object DockerHashResult {
  // See https://docs.docker.com/registry/spec/api/#/content-digests
  val DigestRegex = """([a-zA-Z0-9_+.-]+):([a-zA-Z0-9]+)""".r
  
  def apply(str: String): DockerHashResult = {
    str match {
      case DigestRegex(alg, hash) => new DockerHashResult(alg, hash)
      case _ => throw new IllegalArgumentException(s"Hash value $str does not have the expected 'algorithm:hash' syntax")
    }
  }
}

case class DockerHashResult(hashAlgorithm: String, hashValue: String) {
  val algorithmAndHash = s"$hashAlgorithm:$hashValue"
}
