package wom.values

import java.io.FileNotFoundException
import java.nio.file.NoSuchFileException

import cats.effect.IO
import cats.instances.list._
import cats.syntax.parallel._
import common.validation.IOChecked
import common.validation.IOChecked.IOChecked
import wom.expression.IoFunctionSet
import wom.types._

import scala.util.{Success, Try}

sealed trait WomFile extends WomValue {
  def value: String

  def womFileType: WomFileType

  final def womType: WomType = womFileType

  override def valueString = value

  /**
    * Converts the location using f() recursively converting any files referred to by this file.
    *
    * A caller may need a modification to the path, for example:
    *
    * - Convert a cloud path to a relative path for command line instantiation.
    * - Convert a cloud path to a path on a mounted disk.
    *
    * WomFile such as WomMaybePopulatedFile and WomMaybeListedDirectory may contain references to other WomFile. mapFile
    * will traverse any referenced files, also converting the values within those files too.
    *
    * @param f The function to update the location.
    * @return A new WomFile with the updated location.
    * @see [[wom.values.WomValue.collectAsSeq]]
    * @see [[wom.WomFileMapper.mapWomFiles]]
    */
  def mapFile(f: String => String): WomFile

  /**
    * Converts the location using f() recursively converting any files referred to by this file.
    *
    * A caller may need a modification to the path, for example:
    *
    * - Convert a cloud path to a relative path for command line instantiation.
    * - Convert a cloud path to a path on a mounted disk.
    *
    * WomFile such as WomMaybePopulatedFile and WomMaybeListedDirectory may contain references to other WomFile. mapFile
    * will traverse any referenced files, also converting the values within those files too.
    *
    * Unlike mapFile, mapWomFile is being passed the original WomFile object which allows for more fine tuned mapping
    * based on the type of WomFile for example.
    *
    * @param f The function to update the location.
    * @return A new WomFile with the updated location.
    * @see [[wom.values.WomValue.collectAsSeq]]
    * @see [[wom.WomFileMapper.mapWomFiles]]
    */
  def mapWomFile(f: WomFile => String): WomFile

  /**
    * Converts a file if defined at the partial function.
    *
    * A caller may need a modification to the file, for example:
    *
    * - Fill in the size of the file before the path is localized.
    *
    * WomFile such as WomMaybePopulatedFile and WomMaybeListedDirectory may contain references to other WomFile.
    * mapPartial will traverse any referenced files, also converting the values within those files too.
    *
    * @param f The function to update the location.
    * @return A new WomFile with the updated location.
    * @see [[wom.values.WomValue.collectAsSeq]]
    * @see [[wom.WomFileMapper.mapWomFiles]]
    */
  def collect(f: PartialFunction[WomFile, WomFile]): WomFile

  /**
    * Returns the WomPrimitiveFile instances recursively referenced by this instance.
    *
    * WomMaybeListedDirectory instances return either just the directory as an WomUnlistedDirectory, or if there is a
    * listing then returns the recursive listing. WomMaybePopulatedFile instances return that instance as a
    * WomSingleFile plus any primitives recursively discovered in secondaryFiles.
    *
    * WomPrimitiveFile instances return just the instance.
    */
  def flattenFiles: Seq[WomPrimitiveFile] = {
    this match {
      case womMaybeListedDirectory: WomMaybeListedDirectory =>
        womMaybeListedDirectory.listingOption.getOrElse(Nil).toList match {
          case Nil => womMaybeListedDirectory.valueOption.toList.map(WomUnlistedDirectory)
          case list => list.flatMap(_.flattenFiles)
        }
      case womMaybePopulatedFile: WomMaybePopulatedFile =>
        val primaryFiles: Seq[WomPrimitiveFile] = womMaybePopulatedFile.valueOption.toList.map(WomSingleFile)
        womMaybePopulatedFile.secondaryFiles.foldLeft(primaryFiles) {
          (womFiles, secondaryFile) =>
            womFiles ++ secondaryFile.flattenFiles
        }
      case womPrimitiveFile: WomPrimitiveFile => List(womPrimitiveFile)
    }
  }

  /**
    * If relevant, load the size of the file.
    */
  def withSize(ioFunctionSet: IoFunctionSet): IO[WomFile] = IO.pure(this)
  
  def sizeOption: Option[Long] = None

  protected def recoverFileNotFound[A](fallback: A)(t: Throwable): IO[A] = {
    if (isFileNotFound(t)) IO.pure(fallback) else IO.raiseError(t)
  }

  def isFileNotFound(t: Throwable): Boolean = t match {
    case _: NoSuchFileException | _: FileNotFoundException => true
    case other if other.getCause != null => isFileNotFound(t.getCause)
    case _ => false
  }
}

object WomFile {
  def apply(fileType: WomFileType, value: String) = {
    fileType match {
      case WomUnlistedDirectoryType => WomUnlistedDirectory(value)
      case WomSingleFileType => WomSingleFile(value)
      case WomGlobFileType => WomGlobFile(value)
      case WomMaybeListedDirectoryType => WomMaybeListedDirectory(value)
      case WomMaybePopulatedFileType => WomMaybePopulatedFile(value)
    }
  }
}

sealed trait WomPrimitiveFile extends WomFile with WomPrimitive

/**
  * A directory represented only by a path to a directory.
  *
  * Should not be passed into command line generation. Instead, the execution engine should create a WomListedDirectory
  * locating the files/directories within the `value` and filling in the listing.
  *
  * @param value The location of the directory, possibly in the cloud.
  */
final case class WomUnlistedDirectory(value: String) extends WomPrimitiveFile {
  override val womFileType: WomFileType = WomUnlistedDirectoryType

  override def toWomString = s""""$value""""

  override def add(rhs: WomValue): Try[WomValue] = rhs match {
    case r: WomString => Success(this.copy(value = value + r.value))
    case r: WomOptionalValue => evaluateIfDefined("+", r, add)
    case _ => invalid(s"$value + $rhs")
  }

  override def equals(rhs: WomValue): Try[WomBoolean] = rhs match {
    case r: WomUnlistedDirectory => Success(WomBoolean(this.equals(r)))
    case r: WomString => Success(WomBoolean(value.equals(r.value)))
    case r: WomOptionalValue => evaluateIfDefined("==", r, equals)
    case _ => invalid(s"$value == $rhs")
  }

  override def mapFile(f: String => String): WomUnlistedDirectory = {
    this.copy(value = f(value))
  }

  override def collect(f: PartialFunction[WomFile, WomFile]): WomFile = {
    f.applyOrElse[WomFile, WomFile](this, identity)
  }

  override def mapWomFile(f: WomFile => String) = this.copy(value = f(this))
}

/**
  * A file with no additional files.
  *
  * @param value The location of the file, possibly in the cloud.
  */
final case class WomSingleFile(value: String) extends WomPrimitiveFile {
  override val womFileType: WomFileType = WomSingleFileType

  override def toWomString = s""""$value""""

  override def add(rhs: WomValue): Try[WomValue] = rhs match {
    case r: WomString => Success(this.copy(value = value + r.value))
    case r: WomOptionalValue => evaluateIfDefined("+", r, add)
    case _ => invalid(s"$value + $rhs")
  }

  override def equals(rhs: WomValue): Try[WomBoolean] = rhs match {
    case r: WomSingleFile => Success(WomBoolean(this.equals(r)))
    case r: WomString => Success(WomBoolean(value.equals(r.value)))
    case r: WomOptionalValue => evaluateIfDefined("==", r, equals)
    case _ => invalid(s"$value == $rhs")
  }

  override def mapFile(f: String => String): WomSingleFile = {
    this.copy(value = f(value))
  }

  override def mapWomFile(f: WomFile => String) = {
    this.copy(value = f(this))
  }

  override def collect(f: PartialFunction[WomFile, WomFile]): WomFile = {
    f.applyOrElse[WomFile, WomFile](this, identity)
  }

  override def withSize(ioFunctionSet: IoFunctionSet): IO[WomFile] = {
    IO.fromFuture(IO { ioFunctionSet.size(value)})
      .map(Option.apply)
      .handleErrorWith(recoverFileNotFound(None))
      .map(s => WomMaybePopulatedFile(valueOption = Option(value), sizeOption = s))
  }
}

/**
  * A glob that will be expanded into an array of files from the path in value.
  *
  * Ex:
  * {{{
  *   Array[File] myBams = glob("outdir/\*.bam")
  * }}}
  *
  * @param value The path of the glob within the container.
  */
final case class WomGlobFile(value: String) extends WomPrimitiveFile {
  override val womFileType: WomFileType = WomGlobFileType

  override def toWomString = s"""glob("$value")"""

  override def add(rhs: WomValue): Try[WomValue] = rhs match {
    case r: WomString => Success(this.copy(value + r.value))
    case r: WomOptionalValue => evaluateIfDefined("+", r, add)
    case _ => invalid(s"$value + $rhs")
  }

  override def equals(rhs: WomValue): Try[WomBoolean] = rhs match {
    case r: WomGlobFile => Success(WomBoolean(value.equals(r.value) && womType.equals(r.womType)))
    case r: WomString => Success(WomBoolean(value.toString.equals(r.value.toString)))
    case r: WomOptionalValue => evaluateIfDefined("==", r, equals)
    case _ => invalid(s"$value == $rhs")
  }

  override def mapFile(f: String => String): WomGlobFile = this.copy(value = f(value))

  override def mapWomFile(f: WomFile => String) = this.copy(value = f(this))

  override def collect(f: PartialFunction[WomFile, WomFile]): WomFile = {
    f.applyOrElse[WomFile, WomFile](this, identity)
  }
}


/**
  * A directory possibly with a listing of other files/directories.
  *
  * @param valueOption   The location of the directory, possibly in the cloud.
  * @param listingOption An optional listing of files/directories, either supplied by a user or generated by the engine.
  */
case class WomMaybeListedDirectory(valueOption: Option[String] = None,
                                   listingOption: Option[Seq[WomFile]] = None,
                                   basename: Option[String] = None,
                                   initializeFunction: WomMaybeListedDirectory => IoFunctionSet => IOChecked[WomValue] = { dir => _ => IOChecked.pure(dir) }) extends WomFile {
  override def value: String = {
    valueOption.getOrElse(throw new UnsupportedOperationException(s"value is not available: $this"))
  }

  override val womFileType: WomFileType = WomMaybeListedDirectoryType

  // TODO: WOM: WOMFILE: This isn't even close to a WDL representation (and maybe belongs in WDL?) of this class, but w/o it other areas of the code crash
  override def toWomString = s""""$value""""

  override def mapFile(f: String => String): WomMaybeListedDirectory = {
    this.copy(valueOption = valueOption.map(f), listingOption.map(_.map(_.mapFile(f))))
  }

  override def mapWomFile(f: WomFile => String) = {
    this.copy(valueOption = Option(f(this)), listingOption.map(_.map(_.mapWomFile(f))))
  }

  override def collect(f: PartialFunction[WomFile, WomFile]): WomFile = {
    val copy = this.copy(listingOption = listingOption.map(_.map(_.collect(f))))
    f.applyOrElse[WomFile, WomFile](copy, identity)
  }

  override def withSize(ioFunctionSet: IoFunctionSet): IO[WomFile] = {
    import ioFunctionSet.cs

    listingOption.map({
      _.toList
        .parTraverse[IO, IO.Par, WomFile](_.withSize(ioFunctionSet))
        .map(listingWithSize => this.copy(listingOption = Option(listingWithSize)))
    })
      .getOrElse(IO.pure(this))
  }

  override def initialize(ioFunctionSet: IoFunctionSet): IOChecked[WomValue] = initializeFunction(this)(ioFunctionSet)
  
  override def sizeOption = listingOption.map(_.flatMap(_.sizeOption).sum)
}

object WomMaybeListedDirectory {
  def apply(value: String): WomMaybeListedDirectory = WomMaybeListedDirectory(valueOption = Option(value))
}

/**
  * A file possibly populated with a path plus optional checksum/size/etc.
  *
  * @param valueOption    The location of the file, possibly in the cloud.
  * @param checksumOption An optional checksum of the file contents.
  * @param sizeOption     An optional size of the file contents in bytes.
  * @param formatOption   An optional format description of the file contents.
  * @param contentsOption The optional text contents of the file.
  * @param secondaryFiles Any files associated with this file.
  */
case class WomMaybePopulatedFile(valueOption: Option[String] = None,
                                 checksumOption: Option[String] = None,
                                 override val sizeOption: Option[Long] = None,
                                 formatOption: Option[String] = None,
                                 contentsOption: Option[String] = None,
                                 secondaryFiles: Seq[WomFile] = Vector.empty,
                                 initializeFunction: WomMaybePopulatedFile => IoFunctionSet => IOChecked[WomValue] = { file => _ => IOChecked.pure(file) }
                                ) extends WomFile {
  override def value: String = {
    valueOption.getOrElse(throw new UnsupportedOperationException(s"value is not available: $this"))
  }

  override val womFileType: WomFileType = WomMaybePopulatedFileType

  // TODO: WOM: WOMFILE: This isn't even close to a WDL representation (and maybe belongs in WDL?) of this class, but w/o it other areas of the code crash
  override def toWomString = s""""$value""""

  override def mapFile(f: String => String): WomMaybePopulatedFile = {
    this.copy(valueOption = valueOption.map(f), secondaryFiles = secondaryFiles.map(_.mapFile(f)))
  }

  override def mapWomFile(f: WomFile => String): WomMaybePopulatedFile = {
    this.copy(valueOption = Option(f(this)), secondaryFiles = secondaryFiles.map(_.mapWomFile(f)))
  }

  override def collect(f: PartialFunction[WomFile, WomFile]): WomFile = {
    val copy = this.copy(secondaryFiles = secondaryFiles.map(_.collect(f)))
    f.applyOrElse[WomFile, WomFile](copy, identity)
  }

  override def withSize(ioFunctionSet: IoFunctionSet): IO[WomMaybePopulatedFile] = {
    import ioFunctionSet.cs

    val ioSize: IO[Option[Long]] = (sizeOption, contentsOption) match {
      case (Some(s), _) => IO.pure(Option(s))
      case (None, Some(contents)) => IO.pure(Option(contents.length.toLong))
      case _ => IO.fromFuture(IO { ioFunctionSet.size(value) }).map(Option.apply).handleErrorWith(recoverFileNotFound(sizeOption))
    }
    
    for {
      size <- ioSize
      secondaryFilesWithSize <- secondaryFiles.toList.parTraverse[IO, IO.Par, WomFile](_.withSize(ioFunctionSet))
    } yield this.copy(sizeOption = size, secondaryFiles = secondaryFilesWithSize)
  }

  override def initialize(ioFunctionSet: IoFunctionSet): IOChecked[WomValue] = initializeFunction(this)(ioFunctionSet)
}

object WomMaybePopulatedFile {
  def apply(value: String): WomMaybePopulatedFile = WomMaybePopulatedFile(valueOption = Option(value))
}
