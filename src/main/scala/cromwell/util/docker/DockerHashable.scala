package cromwell.util.docker

sealed trait DockerHashable {
  def getHash: String = ???
}

sealed trait DockerImageId extends DockerHashable

object DockerRegistryImageId {
  val RegistryImageIDLength = 64
}

case class DockerRegistryImageId(id: String) extends DockerImageId

case class DockerHubPartialLayerId(id: String)

// NOTE: Docker Hub returns a list of partial layer IDs back.
case class DockerHubImageId(layerIds: Seq[DockerHubPartialLayerId]) extends DockerImageId

case class DockerFsLayer(blobSum: DockerBlobSum)

case class DockerManifest(fsLayers: Seq[DockerFsLayer]) extends DockerHashable

case class DockerDigestHashable(digest: String) extends DockerHashable
