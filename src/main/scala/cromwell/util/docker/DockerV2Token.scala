package cromwell.util.docker

// https://docs.docker.com/registry/spec/auth/token/#requesting-a-token

/** The properties used to request a token. */
case class DockerV2TokenRequest(realm: String, service: String, scope: Option[String])

/** The returned token. */
case class DockerV2TokenResponse(token: DockerV2Token)
