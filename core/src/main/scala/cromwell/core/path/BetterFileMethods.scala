package cromwell.core.path

import java.io._
import java.net.URI
import java.nio.channels.{AsynchronousFileChannel, FileChannel}
import java.nio.file.attribute._
import java.nio.file.{FileSystem, Files, PathMatcher}
import java.time.Instant
import java.util.zip.Deflater

import scala.collection.JavaConverters._
import scala.io.{BufferedSource, Codec}

/**
  * Implements methods with the same names and signatures as better.files.File.
  *
  * @see [[cromwell.core.path.Path]]
  * @see [[cromwell.core.path.NioPathMethods]]
  * @see [[cromwell.core.path.EvenBetterPathMethods]]
  */
trait BetterFileMethods {
  self: Path =>

  import BetterFileMethods._

  private final def newPath(file: better.files.File): Path = newPath(file.path)

  private final def newPathOrNull(file: better.files.File): Path = Option(file).map(newPath).orNull

  final def toJava: JFile = betterFile.toJava

  final def name: String = betterFile.name

  final def nameOption: Option[String] = betterFile.nameOption

  // NOTE: Check for null before it reaches better files' ==> root: File = [implicit File.apply]path.getRoot
  final def root: Path = Option(betterFile.path.getRoot).map(_ => newPathOrNull(betterFile.root)).orNull

  final def nameWithoutExtension: String = betterFile.nameWithoutExtension

  final def extension: Option[String] = betterFile.extension

  final def extension(includeDot: Boolean = true, includeAll: Boolean = false,
                      toLowerCase: Boolean = true): Option[String] =
    betterFile.extension(includeDot, includeAll, toLowerCase)

  final def hasExtension: Boolean = betterFile.hasExtension

  final def changeExtensionTo(extension: String): Path = newPath(betterFile.changeExtensionTo(extension))

  final def contentType: Option[String] = betterFile.contentType

  final def parent: Path = parentOption.orNull

  final def parentOption: Option[Path] = betterFile.parentOption.map(newPath)

  final def /(child: String): Path = newPath(betterFile./(child))

  final def createChild(child: String, asDirectory: Boolean = false)
                       (implicit attributes: Attributes = Attributes.default,
                        linkOptions: LinkOptions = LinkOptions.default): Path =
    newPath(betterFile.createChild(child, asDirectory)(attributes, linkOptions))

  final def createIfNotExists(asDirectory: Boolean = false, createParents: Boolean = false)
                             (implicit attributes: Attributes = Attributes.default,
                              linkOptions: LinkOptions = LinkOptions.default): Path = {
    newPath(betterFile.createIfNotExists(asDirectory, createParents)(attributes, linkOptions))
  }

  final def exists(implicit linkOptions: LinkOptions = LinkOptions.default): Boolean = betterFile.exists(linkOptions)

  final def notExists(implicit linkOptions: LinkOptions = LinkOptions.default): Boolean =
    betterFile.notExists(linkOptions)

  final def sibling(name: String): Path = newPath(betterFile.sibling(name))

  final def isSiblingOf(sibling: Path): Boolean = betterFile.isSiblingOf(sibling.betterFile)

  final def siblings: Paths = betterFile.siblings.map(newPath)

  final def isChildOf(parent: Path): Boolean = betterFile.isChildOf(parent.betterFile)

  final def contains(file: Path): Boolean = betterFile.contains(file.betterFile)

  final def isParentOf(child: Path): Boolean = contains(child)

  final def bytes: Iterator[Byte] = betterFile.bytes

  final def loadBytes: Array[Byte] = betterFile.loadBytes

  final def byteArray: Array[Byte] = betterFile.byteArray

  final def createDirectory()(implicit attributes: Attributes = Attributes.default): this.type = {
    betterFile.createDirectory()(attributes)
    this
  }

  final def createDirectories()(implicit attributes: Attributes = Attributes.default): this.type = {
    betterFile.createDirectories()(attributes)
    this
  }

  final def chars(implicit codec: Codec): Iterator[Char] = betterFile.chars(codec)

  final def lines(implicit codec: Codec): Traversable[String] = betterFile.lines(codec)

  final def lineIterator(implicit codec: Codec): Iterator[String] = betterFile.lineIterator(codec)

  final def tokens(implicit config: Scanner.Config = Scanner.Config.default, codec: Codec): Traversable[String] =
    betterFile.tokens(config, codec)

  final def contentAsString(implicit codec: Codec): String = betterFile.contentAsString(codec)

  final def `!`(implicit codec: Codec): String = betterFile.`!`(codec)

  final def printLines(lines: Iterator[Any])(implicit openOptions: OpenOptions = OpenOptions.append): this.type = {
    betterFile.printLines(lines)(openOptions)
    this
  }

  final def appendLines(lines: String*)(implicit openOptions: OpenOptions = OpenOptions.append,
                                        codec: Codec): this.type = {
    betterFile.appendLines(lines: _*)(openOptions, codec)
    this
  }

  final def <<(line: String)(implicit openOptions: OpenOptions = OpenOptions.append, codec: Codec): this.type = {
    betterFile.<<(line)(openOptions, codec)
    this
  }

  final def >>:(line: String)(implicit openOptions: OpenOptions = OpenOptions.append, codec: Codec): this.type = {
    betterFile.>>:(line)(openOptions, codec)
    this
  }

  final def appendLine(line: String = "")(implicit openOptions: OpenOptions = OpenOptions.append,
                                          codec: Codec): this.type = {
    betterFile.appendLine(line)(openOptions, codec)
    this
  }

  final def append(text: String)(implicit openOptions: OpenOptions = OpenOptions.append, codec: Codec): this.type = {
    betterFile.append(text)(openOptions, codec)
    this
  }

  final def append(bytes: Array[Byte])(implicit openOptions: OpenOptions): this.type = {
    betterFile.append(bytes)(openOptions)
    this
  }

  final def appendBytes(bytes: Iterator[Byte])(implicit openOptions: OpenOptions = OpenOptions.append): this.type = {
    betterFile.appendBytes(bytes)(openOptions)
    this
  }

  final def write(bytes: Array[Byte])(implicit openOptions: OpenOptions): this.type = {
    betterFile.write(bytes)(openOptions)
    this
  }

  final def writeBytes(bytes: Iterator[Byte])(implicit openOptions: OpenOptions = OpenOptions.append): this.type = {
    betterFile.writeBytes(bytes)(openOptions)
    this
  }

  final def write(text: String)(implicit openOptions: OpenOptions = OpenOptions.default, codec: Codec): this.type = {
    betterFile.write(text)(openOptions, codec)
    this
  }

  final def overwrite(text: String)(implicit openOptions: OpenOptions = OpenOptions.default,
                                    codec: Codec): this.type = {
    betterFile.overwrite(text)(openOptions, codec)
    this
  }

  final def <(text: String)(implicit openOptions: OpenOptions = OpenOptions.default, codec: Codec): this.type = {
    betterFile.<(text)(openOptions, codec)
    this
  }

  final def `>:`(text: String)(implicit openOptions: OpenOptions = OpenOptions.default, codec: Codec): this.type = {
    betterFile.`>:`(text)(openOptions, codec)
    this
  }

  final def newBufferedSource(implicit codec: Codec): BufferedSource = betterFile.newBufferedSource(codec)

  final def bufferedSource(implicit codec: Codec): ManagedResource[BufferedSource] =
    betterFile.bufferedSource(codec)

  final def newRandomAccess(mode: RandomAccessMode = RandomAccessMode.read): RandomAccessFile =
    betterFile.newRandomAccess(mode)

  final def randomAccess(mode: RandomAccessMode = RandomAccessMode.read): ManagedResource[RandomAccessFile] =
    betterFile.randomAccess(mode)

  final def newBufferedReader(implicit codec: Codec): BufferedReader =
    betterFile.newBufferedReader(codec)

  final def bufferedReader(implicit codec: Codec): ManagedResource[BufferedReader] =
    betterFile.bufferedReader(codec)

  final def newBufferedWriter(implicit codec: Codec, openOptions: OpenOptions = OpenOptions.default): BufferedWriter =
    betterFile.newBufferedWriter(codec)

  final def bufferedWriter(implicit codec: Codec,
                           openOptions: OpenOptions = OpenOptions.default): ManagedResource[BufferedWriter] =
    betterFile.bufferedWriter(codec)

  final def newFileReader: FileReader = betterFile.newFileReader

  final def fileReader: ManagedResource[FileReader] = betterFile.fileReader

  final def newFileWriter(append: Boolean = false): FileWriter = betterFile.newFileWriter(append)

  final def fileWriter(append: Boolean = false): ManagedResource[FileWriter] = betterFile.fileWriter(append)

  final def newPrintWriter(autoFlush: Boolean = false)
                          (implicit openOptions: OpenOptions = OpenOptions.default): PrintWriter =
    betterFile.newPrintWriter(autoFlush)

  final def printWriter(autoFlush: Boolean = false)
                       (implicit openOptions: OpenOptions = OpenOptions.default): ManagedResource[PrintWriter] =
    betterFile.printWriter(autoFlush)

  final def newInputStream(implicit openOptions: OpenOptions = OpenOptions.default): InputStream =
    betterFile.newInputStream(openOptions)

  final def inputStream(implicit openOptions: OpenOptions = OpenOptions.default): ManagedResource[InputStream] =
    betterFile.inputStream(openOptions)

  final def newScanner(implicit config: Scanner.Config = Scanner.Config.default): Scanner =
    betterFile.newScanner(config)

  final def scanner(implicit config: Scanner.Config = Scanner.Config.default): ManagedResource[Scanner] =
    betterFile.scanner(config)

  final def newOutputStream(implicit openOptions: OpenOptions = OpenOptions.default): OutputStream =
    betterFile.newOutputStream(openOptions)

  final def outputStream(implicit openOptions: OpenOptions = OpenOptions.default): ManagedResource[OutputStream] =
    betterFile.outputStream(openOptions)

  final def newFileChannel(implicit openOptions: OpenOptions = OpenOptions.default,
                           attributes: Attributes = Attributes.default): FileChannel =
    betterFile.newFileChannel(openOptions, attributes)

  final def fileChannel(implicit openOptions: OpenOptions = OpenOptions.default,
                        attributes: Attributes = Attributes.default): ManagedResource[FileChannel] =
    betterFile.fileChannel(openOptions, attributes)

  final def newAsynchronousFileChannel(implicit openOptions: OpenOptions = OpenOptions.default):
  AsynchronousFileChannel = betterFile.newAsynchronousFileChannel(openOptions)

  final def asynchronousFileChannel(implicit openOptions: OpenOptions = OpenOptions.default):
  ManagedResource[AsynchronousFileChannel] = betterFile.asynchronousFileChannel(openOptions)

  final def digest(algorithmName: String): Array[Byte] = betterFile.digest(algorithmName)

  final def checksum(algorithm: String): String = betterFile.checksum(algorithm)

  final def md5: String = betterFile.md5

  final def sha1: String = betterFile.sha1

  final def sha256: String = betterFile.sha256

  final def sha512: String = betterFile.sha512

  // betterFile.symbolicLink calls Files.readSymbolicLink, but then implicitly converts the java.nio.Path returned to a better.File
  // which calls toAbsolutePath. Consequently, if the path was relative, the current directory is used to make it absolute.
  // This is not the desired behaviour to be able to follow relative symbolic links, so bypass better files method and directly use the java one.
  final def symbolicLink: Option[Path] = {
    if (betterFile.isSymbolicLink) {
      Option(newPath(Files.readSymbolicLink(betterFile.path)))
    } else None
  }

  final def isDirectory: Boolean = betterFile.isDirectory

  final def isRegularFile: Boolean = betterFile.isRegularFile

  final def isSymbolicLink: Boolean = betterFile.isSymbolicLink

  final def isHidden: Boolean = betterFile.isHidden

  final def isLocked(mode: RandomAccessMode, position: Long = 0L, size: Long = Long.MaxValue,
                     isShared: Boolean = false): Boolean = betterFile.isLocked(mode, position, size, isShared)

  final def isReadLocked(position: Long = 0L, size: Long = Long.MaxValue, isShared: Boolean = false): Boolean =
    betterFile.isReadLocked(position, size, isShared)

  final def isWriteLocked(position: Long = 0L, size: Long = Long.MaxValue, isShared: Boolean = false): Boolean =
    betterFile.isWriteLocked(position, size, isShared)

  final def list: Paths = betterFile.list.map(newPath)

  final def children: Paths = betterFile.children.map(newPath)

  final def entries: Paths = betterFile.entries.map(newPath)

  final def listRecursively: Paths = betterFile.listRecursively.map(newPath)

  final def walk(maxDepth: Int = Int.MaxValue): Paths = betterFile.walk(maxDepth).map(newPath)

  final def pathMatcher(syntax: PathMatcherSyntax)(pattern: String): PathMatcher =
    betterFile.pathMatcher(syntax)(pattern)

  final def glob(pattern: String): Paths = betterFile.glob(pattern).map(newPath)

  final def collectChildren(matchFilter: Path => Boolean): Paths = walk().filter(matchFilter(_))

  final def fileSystem: FileSystem = betterFile.fileSystem

  final def uri: URI = betterFile.uri

  final def size: Long = betterFile.size

  final def permissions: Set[PosixFilePermission] = betterFile.permissions

  final def permissionsAsString: String = betterFile.permissionsAsString

  final def setPermissions(permissions: Set[PosixFilePermission]): this.type = {
    betterFile.setPermissions(permissions)
    this
  }

  final def addPermission(permission: PosixFilePermission): this.type = {
    betterFile.addPermission(permission)
    this
  }

  final def removePermission(permission: PosixFilePermission): this.type = {
    betterFile.removePermission(permission)
    this
  }

  // Conflicts with the legacy cromwell.core.path.Obsolete.PathMethodAliases.getFileName(). Uncomment when that's gone.
  //final def apply(permission: PosixFilePermission): Boolean = betterFile.apply(permission)

  final def isOwnerReadable: Boolean = betterFile.isOwnerReadable

  final def isOwnerWritable: Boolean = betterFile.isOwnerWritable

  final def isOwnerExecutable: Boolean = betterFile.isOwnerExecutable

  final def isGroupReadable: Boolean = betterFile.isGroupReadable

  final def isGroupWritable: Boolean = betterFile.isGroupWritable

  final def isGroupExecutable: Boolean = betterFile.isGroupExecutable

  final def isOtherReadable: Boolean = betterFile.isOtherReadable

  final def isOtherWritable: Boolean = betterFile.isOtherWritable

  final def isOtherExecutable: Boolean = betterFile.isOtherExecutable

  final def isReadable: Boolean = betterFile.isReadable

  final def isWriteable: Boolean = betterFile.isWriteable

  final def isExecutable: Boolean = betterFile.isExecutable

  final def attributes(implicit linkOptions: LinkOptions = LinkOptions.default): BasicFileAttributes =
    betterFile.attributes(linkOptions)

  final def posixAttributes(implicit linkOptions: LinkOptions = LinkOptions.default): PosixFileAttributes =
    betterFile.posixAttributes(linkOptions)

  final def dosAttributes(implicit linkOptions: LinkOptions = LinkOptions.default): DosFileAttributes =
    betterFile.dosAttributes(linkOptions)

  final def owner(implicit linkOptions: LinkOptions = LinkOptions.default): UserPrincipal =
    betterFile.owner(linkOptions)

  final def ownerName(implicit linkOptions: LinkOptions = LinkOptions.default): String =
    betterFile.ownerName(linkOptions)

  final def group(implicit linkOptions: LinkOptions = LinkOptions.default): GroupPrincipal =
    betterFile.group(linkOptions)

  final def groupName(implicit linkOptions: LinkOptions = LinkOptions.default): String =
    betterFile.groupName(linkOptions)

  final def setOwner(owner: String): this.type = {
    betterFile.setOwner(owner)
    this
  }

  final def setGroup(group: String): this.type = {
    betterFile.setGroup(group)
    this
  }

  final def touch(time: Instant = Instant.now())
                 (implicit attributes: Attributes = Attributes.default,
                  linkOptions: LinkOptions = LinkOptions.default): this.type = {
    betterFile.touch(time)(attributes, linkOptions)
    this
  }

  final def lastModifiedTime(implicit linkOptions: LinkOptions = LinkOptions.default): Instant =
    betterFile.lastModifiedTime(linkOptions)

  final def delete(swallowIOExceptions: Boolean = false): this.type = {
    betterFile.delete(swallowIOExceptions)
    this
  }

  final def renameTo(newName: String): Path = newPath(betterFile.renameTo(newName))

  final def moveTo(destination: Path, overwrite: Boolean = false): destination.type = {
    betterFile.moveTo(destination.betterFile, overwrite)
    destination
  }

  final def copyTo(destination: Path, overwrite: Boolean = false): destination.type = {
    betterFile.copyTo(destination.betterFile, overwrite)
    destination
  }

  final def symbolicLinkTo(destination: Path)
                          (implicit attributes: Attributes = Attributes.default): destination.type = {
    betterFile.symbolicLinkTo(destination.betterFile)(attributes)
    destination
  }

  final def linkTo(destination: Path, symbolic: Boolean = false)
                  (implicit attributes: Attributes = Attributes.default): destination.type = {
    betterFile.linkTo(destination.betterFile, symbolic)(attributes)
    destination
  }

  final def listRelativePaths(implicit visitOptions: VisitOptions = VisitOptions.default): Paths =
    betterFile.listRelativePaths(visitOptions).map(newPath)

  final def isSamePathAs(that: Path): Boolean = betterFile.isSamePathAs(that.betterFile)

  final def isSameFileAs(that: Path): Boolean = betterFile.isSameFileAs(that.betterFile)

  final def isSameContentAs(that: Path): Boolean = betterFile.isSameContentAs(that.betterFile)

  final def `===`(that: Path): Boolean = betterFile.isSimilarContentAs(that.betterFile)

  final def isSimilarContentAs(that: Path): Boolean = betterFile.isSimilarContentAs(that.betterFile)

  final def =!=(that: Path): Boolean = betterFile.=!=(that.betterFile)

  final def isEmpty: Boolean = betterFile.isEmpty

  final def clear(): this.type = {
    betterFile.clear()
    this
  }

  final def zipTo(destination: Path, compressionLevel: Int = Deflater.DEFAULT_COMPRESSION)
                 (implicit codec: Codec): destination.type = {
    betterFile.zipTo(destination.betterFile, compressionLevel)(codec)
    destination
  }

  final def zip(compressionLevel: Int = Deflater.DEFAULT_COMPRESSION)(implicit codec: Codec): Path =
    newPath(betterFile.zip(compressionLevel)(codec))

  final def unzipTo(destination: Path)(implicit codec: Codec): destination.type = {
    betterFile.unzipTo(destination.betterFile)(codec)
    destination
  }

  final def unzip()(implicit codec: Codec): Path = newPath(betterFile.unzip()(codec))
}

object BetterFileMethods {

  def roots: Iterable[Path] = better.files.File.roots.map(file => DefaultPath(file.path))

  def root: Path = DefaultPath(better.files.File.root.path)

  def home: Path = DefaultPath(better.files.File.home.path)

  def temp: Path = DefaultPath(better.files.File.temp.path)

  def currentWorkingDirectory: Path = DefaultPath(better.files.File.currentWorkingDirectory.path)

  object Cmds {
    def ~ : Path = BetterFileMethods.home

    def pwd: Path = BetterFileMethods.currentWorkingDirectory

    def cwd: Path = pwd

    val `..`: Path => Path = _.parent

    val `.`: Path => Path = identity

    implicit class FileDsl(file: Path) {
      def /(f: Path => Path): Path = f(file)
    }

    def cp(file1: Path, file2: Path): Path = file1.copyTo(file2, overwrite = true)

    def mv(file1: Path, file2: Path): Path = file1.moveTo(file2, overwrite = true)

    def rm(file: Path): Path = file.delete(swallowIOExceptions = true)

    def del(file: Path): Path = rm(file)

    def ln(file1: Path, file2: Path): Path = file1 linkTo file2

    def ln_s(file1: Path, file2: Path): Path = file1 symbolicLinkTo file2

    def cat(files: Path*): Seq[Iterator[Byte]] = files.map(_.bytes)

    def ls(file: Path): Paths = file.list

    def dir(file: Path): Paths = ls(file)

    def ls_r(file: Path): Paths = file.listRecursively

    def touch(file: Path): Path = file.touch()

    def mkdir(file: Path): Path = file.createDirectory()

    def md5(file: Path): String = file.md5

    def sha1(file: Path): String = file.sha1

    def sha256(file: Path): String = file.sha256

    def sha512(file: Path): String = file.sha512

    def mkdirs(file: Path): Path = file.createDirectories()

    def chown(owner: String, file: Path): Path = file.setOwner(owner)

    def chgrp(group: String, file: Path): Path = file.setGroup(group)

    def chmod(permissions: String, file: Path): Path = file.setPermissions(PosixFilePermissions.fromString(permissions).asScala.toSet)

    def chmod_+(permission: PosixFilePermission, file: Path): Path = file.addPermission(permission)

    def chmod_-(permission: PosixFilePermission, file: Path): Path = file.removePermission(permission)

    def stat(file: Path): PosixFileAttributes = file.posixAttributes

    def unzip(zipFile: Path)(destination: Path)(implicit codec: Codec): destination.type =
      zipFile.unzipTo(destination)(codec)

    def zip(files: Path*)(destination: Path, compressionLevel: Int = Deflater.DEFAULT_COMPRESSION)
           (implicit codec: Codec): destination.type = {
      better.files.Cmds.zip(files.map(_.betterFile): _*)(destination.betterFile, compressionLevel)
      destination
    }
  }

  type PathMatcherSyntax = better.files.File.PathMatcherSyntax
  val PathMatcherSyntax = better.files.File.PathMatcherSyntax

  type RandomAccessMode = better.files.File.RandomAccessMode
  val RandomAccessMode = better.files.File.RandomAccessMode

  type Attributes = better.files.File.Attributes
  val Attributes = better.files.File.Attributes

  type CopyOptions = better.files.File.CopyOptions
  val CopyOptions = better.files.File.CopyOptions

  type Events = better.files.File.Events
  val Events = better.files.File.Events

  type OpenOptions = better.files.File.OpenOptions
  val OpenOptions = better.files.File.OpenOptions

  type LinkOptions = better.files.File.LinkOptions
  val LinkOptions = better.files.File.LinkOptions

  type VisitOptions = better.files.File.VisitOptions
  val VisitOptions = better.files.File.VisitOptions

  type Order = better.files.File.Order
  val Order = better.files.File.Order

  type Scanner = better.files.Scanner
  val Scanner = better.files.Scanner

  type JFile = java.io.File
  type Paths = Iterator[Path]

  type Closable = better.files.Closeable
  type ManagedResource[A <: Closable] = better.files.ManagedResource[A]
}
