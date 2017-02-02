package cromwell.core.callcaching.docker

case class DockerHashRequest(dockerImageID: DockerImageIdentifierWithoutHash, credentials: Seq[Any] = Seq.empty)
