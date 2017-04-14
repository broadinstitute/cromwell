package cromwell.docker

case class DockerHashRequest(dockerImageID: DockerImageIdentifierWithoutHash, credentials: List[Any] = List.empty)
