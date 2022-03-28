package cromwell.docker

import cromwell.core.DockerCredentialUsernameAndPassword

case class DockerInfoRequest(dockerImageID: DockerImageIdentifier, credentials: List[Any] = List.empty) {
  def credentialDetails: List[String] = credentials map {
    case DockerCredentialUsernameAndPassword(username, _) => s"DockerCredentials($username, ...)"
    case other => s"${other.getClass.getSimpleName}(...)"
  }
}
