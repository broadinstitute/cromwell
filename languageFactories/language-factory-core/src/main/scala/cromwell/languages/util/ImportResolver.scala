package cromwell.languages.util

import java.net.{URI, URL}
import java.nio.file.Paths

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
import cromwell.core.path.Path
import java.nio.file.{Path => NioPath}

import wom.core.WorkflowSource

import scala.concurrent.duration._
import scala.concurrent.Await
import scala.util.{Failure, Success, Try}

object ImportResolver {

  case class ImportResolutionRequest(toResolve: String, currentResolvers: List[ImportResolver])
  case class ResolvedImportBundle(source: WorkflowSource, newResolvers: List[ImportResolver])

  sealed trait ImportResolver {
    def name: String
    protected def innerResolver(path: String, currentResolvers: List[ImportResolver]): Checked[ResolvedImportBundle]
    def resolver: CheckedAtoB[ImportResolutionRequest, ResolvedImportBundle] = CheckedAtoB.fromCheck { request =>
      innerResolver(request.toResolve, request.currentResolvers).contextualizeErrors(s"resolve '${request.toResolve}' using resolver: '$name'")
    }
  }


  object DirectoryResolver {
    def apply(directory: Path, allowEscapingDirectory: Boolean): DirectoryResolver = {
      val dontEscapeFrom = if (allowEscapingDirectory) None else Option(directory.toJava.getCanonicalPath)
      DirectoryResolver(directory, dontEscapeFrom)
    }
  }

  case class DirectoryResolver(directory: Path, dontEscapeFrom: Option[String] = None) extends ImportResolver {
    lazy val absolutePathToDirectory = directory.toJava.getCanonicalPath

    override def innerResolver(path: String, currentResolvers: List[ImportResolver]): Checked[ResolvedImportBundle] = {

      def updatedResolverSet(oldRootDirectory: Path, newRootDirectory: Path, current: List[ImportResolver]): List[ImportResolver] = {
        current map {
          case d if d == this => DirectoryResolver(newRootDirectory, dontEscapeFrom)
          case other => other
        }
      }

      def fetchContentFromAbsolutePath(absolutePathToFile: NioPath): ErrorOr[String] = {
        checkLocation(absolutePathToFile, path) flatMap { _ =>
          val file = File(absolutePathToFile)
          if (file.exists) {
            File(absolutePathToFile).contentAsString.validNel
          } else {
            s"Import file not found: $path".invalidNel
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

    override def name: String = s"relative to directory $absolutePathToDirectory (without escaping $dontEscapeFrom)"
  }

  def zippedImportResolver(zippedImports: Array[Byte]): ErrorOr[ImportResolver] = {
    LanguageFactoryUtil.validateImportsDirectory(zippedImports) map { dir =>
      DirectoryResolver(dir, Option(dir.toJava.getCanonicalPath))
    }
  }

  case class HttpResolver(relativeTo: Option[String] = None, headers: Map[String, String] = Map.empty) extends ImportResolver {
    import HttpResolver._

    override def name: String = {
      s"http importer${relativeTo map { r => s" (relative to $r)" } getOrElse ""}"
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
        else s"Cannot import '$str' relative to nothing".invalidNelCheck
    }

    override def innerResolver(str: String, currentResolvers: List[ImportResolver]): Checked[ResolvedImportBundle] = {
      pathToLookup(str) flatMap { toLookup =>
        (Try {
          implicit val sttpBackend = HttpResolver.sttpBackend
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
  }

  object HttpResolver {

    lazy val sttpBackend = AsyncHttpClientCatsBackend[IO]()

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
