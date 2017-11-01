package centaur.test

import com.google.cloud.storage.{Blob, Storage}
import com.google.cloud.storage.Storage.BlobListOption

import scala.collection.JavaConverters._

import scala.language.implicitConversions

sealed trait CheckFiles {
  def countObjectsAtPath: String => Int
}

case class JesCheckFiles(storage: Storage) extends CheckFiles {
  import GCS._

  def countObjectsAtPath: String => Int =
    storage.parsePath andThen storage.countObjectsAtPath
}


case class LocalCheckFiles() extends CheckFiles {
  def countObjectsAtPath: String => Int = { s =>
    val d = new java.io.File(s)
    if (d.exists && d.isDirectory)
      d.listFiles.length
    else if (d.exists && d.isFile)
      1
    else
      0
  }
}

case class GCSPath(bucket: String, directory: String)

object GCS {

  implicit def gcsOps(s: Storage): GCSOps = GCSOps(s)
}

case class GCSOps(storage: Storage) {

  def parsePath: String => GCSPath =
    {fullPath =>
      val bucketAndDashes = fullPath.drop(5).split('/')
      val bucket = bucketAndDashes.head
      val directory = bucketAndDashes.tail.mkString("/")

      GCSPath(bucket, directory)
    }

  private def listObjectsAtPath: GCSPath => Iterable[Blob] =
    g =>
      storage.list(g.bucket, BlobListOption.prefix(g.directory)).iterateAll.asScala

  def countObjectsAtPath: GCSPath => Int =
      listObjectsAtPath(_).size
}
