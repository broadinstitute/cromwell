package cloud.nio.spi

/**
  * Contains a listing of paths/keys under some host/bucket. The keys must be absolute paths but should not begin with a
  * slash.
  */
case class CloudNioFileList(paths: Iterable[String], markerOption: Option[String])
