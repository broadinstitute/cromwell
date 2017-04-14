package cromwell.docker.local

/**
  * Used by the docker cli for lookups or pulls.
  *
  * The format of `repository` is based on the output of `docker images`. It is also valid as input for the cli
  * commands.
  *
  * @param repository The repository as listed by `docker images`. Includes the image host and name, but not the tag.
  * @param tag The docker image tag.
  */
case class DockerCliKey(repository: String, tag: String) {
  val fullName = s"$repository:$tag"
}
