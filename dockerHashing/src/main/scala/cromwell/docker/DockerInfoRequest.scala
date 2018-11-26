package cromwell.docker

case class DockerInfoRequest(dockerImageID: DockerImageIdentifier, credentials: List[Any] = List.empty)
