package cromwell.util.docker

/**
  * A set of endpoints belonging to a registry.
  *
  * For most registries, such as gcr.io, quay.io, etc., the namespace, the registry v1 hostname, and the registry v2
  * hostname are all the same.
  *
  * Docker Hub is an example registry that uses a different namespace, a different registry v1 hostname, and a
  * different registry v2 hostname.
  *
  * - https://github.com/docker/docker/blob/v1.9.1/registry/config.go#L24-L25
  * - https://github.com/docker/docker/blob/v1.9.1/registry/config_unix.go#L6-L10
  */
object DockerRegistry {
  /**
    * Creates a registry where the namespace, the registry v1 hostname, and the registry v2 hostname are all the same.
    *
    * @param host The host to use as the namespace and registry endpoints.
    * @param login The login information for the registry.
    * @return The DockerRegistry.
    */
  def apply(host: String, login: DockerLoginProvider): DockerRegistry = DockerRegistry(host, host, host, login)
}

/**
  * Creates a registry.
  *
  * @param namespace The namespace of the registry.
  * @param v1Hostname The host for contacting the V1 API.
  * @param v2Hostname The host for contacting the V2 API.
  */
case class DockerRegistry(namespace: String, v1Hostname: String, v2Hostname: String, loginProvider: DockerLoginProvider)
