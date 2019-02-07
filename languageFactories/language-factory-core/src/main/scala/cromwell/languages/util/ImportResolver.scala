package cromwell.languages.util

import java.net.{URI, URL}
import java.nio.file.{FileSystem, FileSystems, Files, Paths, Path => NioPath}

import com.google.common.jimfs.Configuration
import com.google.common.jimfs.Jimfs
import better.files.File
import cats.data.NonEmptyList
import cats.effect.IO
import cats.syntax.either._
import cats.syntax.validated._
import com.softwaremill.sttp._
import com.softwaremill.sttp.asynchttpclient.cats.AsyncHttpClientCatsBackend
import common.Checked
import common.transforms.CheckedAtoB
import common.validation.ErrorOr._
import common.validation.Checked._
import common.validation.Validation._
import cromwell.core.path.{DefaultPathBuilder, Path}
import cromwell.core.WorkflowId
import wom.core.WorkflowSource

import scala.concurrent.duration._
import scala.concurrent.Await
import scala.util.{Failure, Success, Try}

object ImportResolver {

  case class ImportResolutionRequest(toResolve: String, currentResolvers: List[ImportResolver])
  case class ResolvedImportBundle(source: WorkflowSource, newResolvers: List[ImportResolver])

  trait ImportResolver {
    def name: String
    protected def innerResolver(path: String, currentResolvers: List[ImportResolver]): Checked[ResolvedImportBundle]
    def resolver: CheckedAtoB[ImportResolutionRequest, ResolvedImportBundle] = CheckedAtoB.fromCheck { request =>
      innerResolver(request.toResolve, request.currentResolvers).contextualizeErrors(s"resolve '${request.toResolve}' using resolver: '$name'")
    }
    def close(): Try[Unit]
  }

  object DirectoryResolver {
    def apply(directory: Path, allowEscapingDirectory: Boolean, customName: Option[String]): DirectoryResolver = {
      val dontEscapeFrom = if (allowEscapingDirectory) None else Option(directory.toJava.getCanonicalPath)
      DirectoryResolver(directory, dontEscapeFrom, customName)
    }

    def localFilesystemResolvers(baseWdl: Option[Path]) = List(
      DirectoryResolver(
        DefaultPathBuilder.build(Paths.get(".")),
        allowEscapingDirectory = true,
        customName = None
      ),
      DirectoryResolver(
        DefaultPathBuilder.build(Paths.get("/")),
        allowEscapingDirectory = false,
        customName = Some("entire local filesystem (relative to '/')")
      )
    ) ++ baseWdl.toList.map { rt =>
      DirectoryResolver(
        DefaultPathBuilder.build(Paths.get(rt.toAbsolutePath.toFile.getParent)),
        allowEscapingDirectory = true,
        customName = None
      )
    }
  }

  case class DirectoryResolver(directory: Path,
                               dontEscapeFrom: Option[String] = None,
                               customName: Option[String]) extends ImportResolver {
    lazy val absolutePathToDirectory: String = directory.toJava.getCanonicalPath

    override def innerResolver(path: String, currentResolvers: List[ImportResolver]): Checked[ResolvedImportBundle] = {

      def updatedResolverSet(oldRootDirectory: Path, newRootDirectory: Path, current: List[ImportResolver]): List[ImportResolver] = {
        current map {
          case d if d == this => DirectoryResolver(newRootDirectory, dontEscapeFrom, customName)
          case other => other
        }
      }

      def fetchContentFromAbsolutePath(absolutePathToFile: NioPath): ErrorOr[String] = {
        checkLocation(absolutePathToFile, path) flatMap { _ =>
          val file = File(absolutePathToFile)
          if (file.exists) {
            File(absolutePathToFile).contentAsString.validNel
          } else {
            s"File not found: $path".invalidNel
          }
        }
      }

      val errorOr = for {
        resolvedPath <- resolvePath(path)
        absolutePathToFile <- makeAbsolute(resolvedPath)
        fileContents <- fetchContentFromAbsolutePath(absolutePathToFile)
        updatedResolvers = updatedResolverSet(directory, resolvedPath.parent, currentResolvers)
      } yield ResolvedImportBundle(fileContents, updatedResolvers)

      errorOr.toEither
    }

    private def resolvePath(path: String): ErrorOr[Path] = Try(directory.resolve(path)).toErrorOr
    private def makeAbsolute(resolvedPath: Path): ErrorOr[NioPath] = Try(Paths.get(resolvedPath.toFile.getCanonicalPath)).toErrorOr
    private def checkLocation(absoluteNioPath: NioPath, reportedPathIfBad: String): ErrorOr[Unit] =
      if (dontEscapeFrom.forall(absoluteNioPath.startsWith))
        ().validNel
      else
        s"$reportedPathIfBad is not an allowed import path".invalidNel

    def resolveAndMakeAbsolute(path: String): ErrorOr[NioPath] = for {
      resolved <- resolvePath(path)
      abs <- makeAbsolute(resolved)
      _ <- checkLocation(abs, path)
    } yield abs

    override lazy val name: String = (customName, dontEscapeFrom) match {
      case (Some(custom), _) => custom
      case (None, Some(dontEscapePath)) =>
        val dontEscapeFromDirname = Paths.get(dontEscapePath).getFileName.toString
        val shortPathToDirectory = absolutePathToDirectory.stripPrefix(dontEscapePath)

        val relativePathToDontEscapeFrom = s"[...]/$dontEscapeFromDirname"
        val relativePathToDirectory = s"$relativePathToDontEscapeFrom$shortPathToDirectory"

        s"relative to directory $relativePathToDirectory (without escaping $relativePathToDontEscapeFrom)"
      case (None, None) =>
        val shortPathToDirectory = Paths.get(absolutePathToDirectory).toFile.getCanonicalFile.toPath.getFileName.toString
        s"relative to directory [...]/$shortPathToDirectory (escaping allowed)"
    }

    override def close(): Try[Unit] = Success(())
  }

  // !!! Warning: does not yet work for nested directories !!!
  case class ZipResolver(zipContents: Array[Byte], workflowId: WorkflowId) extends ImportResolver {

    // Based on https://github.com/google/jimfs/issues/47#issuecomment-311360741
    lazy val zipMemoryMappedFilesystem: FileSystem = {

      // Write the zip bytes to a specific path on a fresh Google JimFS in-memory filesystem and return the path
      lazy val zipPath = {
        // Use workflow ID and `lazy` because name must be globally unique
        val zipFilesystem = Jimfs.newFileSystem("jimfs-workflow-" + workflowId.toString, Configuration.unix())

        val zipPath = zipFilesystem.getPath("/imports.zip")
        Files.write(zipPath, zipContents)

        zipPath
      }

      val env = new java.util.HashMap[String, String]()
      env.put("create", "false") // `true` means "create the zip if it doesn't exist" - would mean we have a bug
      env.put("encoding", "UTF-8")

      // Map the zip into memory using the Java zip filesystem provider
      FileSystems.newFileSystem(new java.net.URI("jar:" + zipPath.toUri), env)
    }

    override def name: String = s"Zip Resolver for zip of length ${zipContents.length}"

    override protected def innerResolver(path: String, currentResolvers: List[ImportResolver]): Checked[ResolvedImportBundle] = {
      val pathOnFilesystem: NioPath = zipMemoryMappedFilesystem.getPath(path)
      val contents = File(pathOnFilesystem).contentAsString

      ResolvedImportBundle(contents, List()).validNelCheck
    }

    override def close(): Try[Unit] = Try {
      // Not sure if this is necessary
      zipMemoryMappedFilesystem.close()
    }
  }

  // This is where we commit the original sin of turning into an in-memory Array[Byte] into a stateful global eternal directory on the filesystem
  def zippedImportResolver(zippedImports: Array[Byte], workflowId: WorkflowId): ErrorOr[ZipResolver] = {
    Try(ZipResolver(zippedImports, workflowId)).toErrorOr
  }

  case class HttpResolver(relativeTo: Option[String] = None, headers: Map[String, String] = Map.empty) extends ImportResolver {
    import HttpResolver._

    override def name: String = relativeTo match {
      case Some(relativeToPath) => s"http importer (relative to $relativeToPath)"
      case None => "http importer (no 'relative-to' origin)"

    }

    def newResolverList(newRoot: String): List[ImportResolver] = {
      val rootWithoutFilename = newRoot.split('/').init.mkString("", "/", "/")
      List(
        HttpResolver(relativeTo = Some(canonicalize(rootWithoutFilename)), headers)
      )
    }

    def pathToLookup(str: String): Checked[String] = relativeTo match {
      case Some(relativeToValue) =>
        if (str.startsWith("http")) canonicalize(str).validNelCheck
        else if (str.startsWith("/")) canonicalize(s"${root(relativeToValue).stripSuffix("/")}/$str").validNelCheck
        else canonicalize(s"${relativeToValue.stripSuffix("/")}/$str").validNelCheck
      case None =>
        if (str.startsWith("http")) canonicalize(str).validNelCheck
        else "Relative path".invalidNelCheck
    }

    override def innerResolver(str: String, currentResolvers: List[ImportResolver]): Checked[ResolvedImportBundle] = {
      pathToLookup(str) flatMap { toLookup =>
        (Try {
          implicit val sttpBackend = HttpResolver.sttpBackend()
          val responseIO: IO[Response[String]] = sttp.get(uri"$toLookup").headers(headers).send()

          // temporary situation to get functionality working before
          // starting in on async-ifying the entire WdlNamespace flow
          val result: Checked[String] = Await.result(responseIO.unsafeToFuture, 15.seconds).body.leftMap(NonEmptyList(_, List.empty))

          result map {
            ResolvedImportBundle(_, newResolverList(toLookup))
          }
        } match {
          case Success(result) => result
          case Failure(e) => s"HTTP resolver with headers had an unexpected error (${e.getMessage})".invalidNelCheck
        }).contextualizeErrors(s"download $toLookup")
      }
    }

    override def close(): Try[Unit] = Success(())
  }

  object HttpResolver {

    import common.util.IntrospectableLazy
    import common.util.IntrospectableLazy._

    val sttpBackend: IntrospectableLazy[SttpBackend[IO, Nothing]] = lazily {
      AsyncHttpClientCatsBackend[IO]()
    }

    def closeBackendIfNecessary() = if (sttpBackend.exists) sttpBackend.close()

    def canonicalize(str: String): String = {
      val url = new URL(str)
      val uri = url.toURI.normalize
      uri.toString
    }

    def root(str: String): String = {
      val url = new URL(str)
      val uri = new URI(url.getProtocol, url.getUserInfo, url.getHost, url.getPort, null, null, null).normalize()
      uri.toString
    }
  }

}
