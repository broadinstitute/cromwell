//package cromwell.filesystems.demo.dos
//
//import com.google.api.client.testing.http.MockHttpTransport
//import com.google.api.client.testing.json.MockJsonFactory
//import com.google.api.services.storage.Storage
//import com.google.cloud.storage.StorageOptions
//import com.google.cloud.storage.contrib.nio.{CloudStorageConfiguration, CloudStoragePath}
//import cromwell.core.path.{NioPath, Path, PathBuilder}
//import cromwell.filesystems.demo.dos.DemoDosPathBuilder._
//import cromwell.filesystems.gcs.{GcsPath, GcsPathBuilder}
//
//import scala.util.{Failure, Success, Try}
//
//class DemoDosPathBuilder extends PathBuilder {
//  private val mockStorage: Storage =
//    new Storage.Builder(new MockHttpTransport, new MockJsonFactory, null)
//      .setApplicationName("cromwell-demo-dos")
//      .build()
//
//  private val gcsPathBuilder: GcsPathBuilder = {
//    new GcsPathBuilder(mockStorage, CloudStorageConfiguration.DEFAULT, StorageOptions.getDefaultInstance)
//  }
//
//  override def name: String = "Demo Dos"
//
//  override def build(pathAsString: String): Try[Path] = {
//    if (DemoDosPathBuilder.accepts(pathAsString)) {
//      gcsPathBuilder.build(pathAsString.swapPrefix("dos", "gs")).transform(
//        gcsPath => Success(DemoDosPath(gcsPath, this)),
//        gcsException => {
//          val dosMessage = gcsException.getMessage
//            .replaceAll("gs://", "dos://")
//            .replaceAll("GCS", "DOS")
//          val dosException = new RuntimeException(dosMessage, gcsException)
//          Failure(dosException)
//        }
//      )
//    } else {
//      Failure(new IllegalArgumentException(s"$pathAsString does not have a dos scheme"))
//    }
//  }
//}
//
//object DemoDosPathBuilder {
//
//  implicit class EnhancedString(val string: String) extends AnyVal {
//    def swapPrefix(from: String, to: String) = {
//      if (string.startsWith(from)) {
//        to + string.stripPrefix(from)
//      } else {
//        string
//      }
//    }
//  }
//
//  def accepts(path: String): Boolean = path.startsWith("dos://")
//}
//
//case class DemoDosPath(gcsPath: GcsPath, demoDosPathBuilder: DemoDosPathBuilder)
//  extends Path {
//  override protected def nioPath: NioPath = gcsPath.nioPath
//
//  override protected def newPath(nioPath: NioPath): Path = {
//    nioPath match {
//      case cloudStoragePath: CloudStoragePath =>
//        val host = cloudStoragePath.bucket().stripSuffix("/")
//        val path = cloudStoragePath.toString.stripPrefix("/")
//        demoDosPathBuilder.build(s"dos://$host/$path").get
//      case _ =>
//        throw new RuntimeException(
//          s"Internal path was not a cloud storage path: $nioPath")
//    }
//  }
//
//  override def pathAsString: String = {
//    gcsPath.pathAsString.swapPrefix("gs", "dos")
//  }
//
//  override def pathWithoutScheme: String = gcsPath.pathWithoutScheme
//}
