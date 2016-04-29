package cromwell.util.docker

// https://docs.docker.com/registry/spec/auth/token/#requesting-a-token

/** The properties used to request a token. */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class DockerV2TokenRequest(realm: String, service: String, scope: Option[String])

/** The returned token. */
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class DockerV2TokenResponse(token: DockerV2Token)
