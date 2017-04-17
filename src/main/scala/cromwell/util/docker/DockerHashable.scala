package cromwell.util.docker

/** Anything that can provide a unique hash of the docker image. */
sealed trait DockerHashable {
  /**
    * TBD: hashing implementation for each class, that meets specs for
    * what will be stored in the databse for job avoidance.
    */
  def getHash: String = ???
}

/** Returned by Docker when asked for a V1 image tag, as opposed to a V2 manifest. */
sealed trait DockerImageId extends DockerHashable

// Companion
object DockerRegistryImageId {
  /**
    * Length of the id in DockerRegistryImageId
    */
  val RegistryImageIDLength = 64
}

/**
  * Returned by the deprecated Docker Registry API V1 when asked for an image tag.
  * https://docs.docker.com/v1.6/reference/api/registry_api/#get-image-id-for-a-particular-tag
  */
case class DockerRegistryImageId(id: String) extends DockerImageId

/** The "id" part of the DockerHubImageId object returned by Docker Hub. */
case class DockerHubPartialLayerId(id: String)

/**
  * The list of layer ids returned by Docker Hub when asked for an image tag.
  * NOTE: Docker Hub returns a list of layer IDs truncated to 8 characters, not the full hash of each layer.
  */
case class DockerHubImageId(layerIds: Seq[DockerHubPartialLayerId]) extends DockerImageId

/** Returned as part of the V2 Manifest */
case class DockerFsLayer(blobSum: DockerBlobSum)

/** Returned by docker when asked for a docker V2 manifest, as opposed to an image. */
case class DockerManifest(fsLayers: Seq[DockerFsLayer]) extends DockerHashable

/**
  * Docker supports downloading images via a digest identifier.
  * NOTE: We don't need to check the registry for these hashes, and can just pull them from the user specified
  * identifier.
  */
case class DockerDigestHashable(digest: String) extends DockerHashable
