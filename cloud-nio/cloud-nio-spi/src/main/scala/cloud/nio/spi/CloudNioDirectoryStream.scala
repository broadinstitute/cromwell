package cloud.nio.spi

import java.nio.file.{DirectoryStream, Path}

import scala.jdk.CollectionConverters._

class CloudNioDirectoryStream(
  fileProvider: CloudNioFileProvider,
  retry: CloudNioRetry,
  prefix: CloudNioPath,
  filter: DirectoryStream.Filter[_ >: Path]
) extends DirectoryStream[Path] {

  override def iterator(): java.util.Iterator[Path] = pathStream().filterNot(_ == prefix).iterator.asJava

  private[this] def pathStream(markerOption: Option[String] = None): LazyList[Path] =
    listNext(markerOption) match {
      case CloudNioFileList(keys, Some(marker)) =>
        keys.to(LazyList).map(toPath) ++ pathStream(Option(marker))
      case CloudNioFileList(keys, None) =>
        keys.to(LazyList).map(toPath)
    }

  private[this] def toPath(key: String): Path =
    prefix.getFileSystem.getPath("/" + key)

  private[this] def listNext(markerOption: Option[String]): CloudNioFileList =
    retry.from(() => fileProvider.listObjects(prefix.cloudHost, prefix.cloudPath, markerOption))

  override def close(): Unit = {}

}
