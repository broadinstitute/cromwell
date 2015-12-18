package cromwell.util.docker

import scala.util.Try

/** Anything that can provide a unique hash of the docker image. */
sealed trait DockerHashable {
  /** A descriptor for the hash */
  val hashCollectionType: String

  /**
    * Strings that uniquely identify this docker image, that change if there are any modifications to the image.
    * The string may be hashes themselves, optionally prefixed with a type. The only requirement is that they
    * change when the image type changes.
    */
  val dockerHashes: Seq[Try[DockerHash]]

  /** A unique hash that changes when the docker image changes, prefixed with the type of the generated hash. */
  lazy val dockerHash: Try[DockerHash] = DockerHash.fromTries(hashCollectionType, dockerHashes)
}

/** Returned by Docker when asked for a V1 image tag, as opposed to a V2 manifest. */
sealed trait DockerImageId extends DockerHashable

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
case class DockerRegistryImageId(id: String) extends DockerImageId {
  override val hashCollectionType = "imageId"
  override val dockerHashes = Seq(DockerHash.fromHash("sha256", id))
}

/** The "id" part of the DockerHubImageId object returned by Docker Hub. */
case class DockerHubPartialLayerId(id: String)

/**
  * The list of layer ids returned by Docker Hub when asked for an image tag.
  * NOTE: Docker Hub returns a list of layer IDs truncated to 8 characters, not the full hash of each layer.
  */
case class DockerHubImageId(layerIds: Seq[DockerHubPartialLayerId]) extends DockerImageId {
  override val hashCollectionType = "layerIds"
  override val dockerHashes = layerIds.map(layerId => DockerHash.fromHash("sha256Part", layerId.id))
}

/** Returned as part of the V2 Manifest */
case class DockerFsLayer(blobSum: DockerBlobSum)

/** Returned by docker when asked for a docker V2 manifest, as opposed to an image. */
case class DockerManifest(fsLayers: Seq[DockerFsLayer]) extends DockerHashable {
  override val hashCollectionType = "layerBlobs"
  override val dockerHashes = fsLayers.map(fsLayer => DockerHash.fromDigest(fsLayer.blobSum))
}

/**
  * Docker supports downloading images via a digest identifier.
  * NOTE: We don't need to check the registry for these hashes, and can just pull them from the user specified
  * identifier.
  */
case class DockerDigestHashable(digest: String) extends DockerHashable {
  override val hashCollectionType = "digest"
  override val dockerHashes = Seq(DockerHash.fromDigest(digest))
}
