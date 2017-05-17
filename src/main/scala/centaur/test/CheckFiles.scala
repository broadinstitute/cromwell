package centaur.test

import com.google.cloud.storage.{Storage, Blob, BlobId}
import com.google.cloud.storage.Storage.BlobListOption
import scala.collection.JavaConverters._
import java.io.File

sealed trait CheckFiles {
  def checkDirectorySize: String => Int 
  def deleteExistingFiles: String => Unit
}

case class JesCheckFiles(storage: Storage) extends CheckFiles {
  import GCS._
  import storage._

  def checkDirectorySize: String => Int = 
    storage.parsePath andThen storage.listDirectoryCount

  def deleteExistingFiles: String => Unit = 
    storage.deleteDirectory
}


case class LocalCheckFiles() extends CheckFiles {
  def checkDirectorySize: String => Int = {s =>
    val d = new java.io.File(s)
    if (d.exists && d.isDirectory)
      d.listFiles.size
    else 
      0
  }

  def deleteExistingFiles = 
    {directory: String =>
      val dir = new java.io.File(directory)
      def deleteDir(d: File) {
        if (d.exists && d.isDirectory) {
          d.listFiles.toList.foreach { f =>

          //directories must be empty to be deleted
          if (f.isDirectory) 
            deleteDir(f)

          f.delete
          }
        }
      }

      deleteDir(dir)
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

  def listDirectory: GCSPath => Iterable[Blob] =
    g => 
      storage.list(g.bucket, BlobListOption.prefix(g.directory)).iterateAll.asScala

  def listDirectoryCount: GCSPath => Int =
      listDirectory(_).size

  def deleteDirectory: String => Unit = 
    dir => {
      val gcsPath = parsePath(dir)
      val blobs = listDirectory(gcsPath)
      if (blobs.size > 0) {
        val blobIds = blobs.map(_.getBlobId).asJava
        storage.delete(blobIds)
      }
    }
}
