package cloud.nio.spi

import java.nio.file.{DirectoryStream, Path}

import scala.collection.JavaConverters._

class CloudNioDirectoryStream(
  fileProvider: CloudNioFileProvider,
  retry: CloudNioRetry,
  prefix: CloudNioPath,
  filter: DirectoryStream.Filter[_ >: Path]
) extends DirectoryStream[Path] {

  override def iterator(): java.util.Iterator[Path] = pathStream().filterNot(_ == prefix).toIterator.asJava

  private[this] def pathStream(markerOption: Option[String] = None): Stream[Path] = {
    listNext(markerOption) match {
      case CloudNioFileList(keys, Some(marker)) =>
        keys.toStream.map(toPath) ++ pathStream(Option(marker))
      case CloudNioFileList(keys, None) =>
        keys.toStream.map(toPath)
    }
  }

  private[this] def toPath(key: String): Path = {
    prefix.getFileSystem.getPath("/" + key)
  }

  private[this] def listNext(markerOption: Option[String]): CloudNioFileList = {
    retry.from(
      () => fileProvider.listObjects(prefix.cloudHost, prefix.cloudPath, markerOption)
    )
  }

  override def close(): Unit = {}

}
