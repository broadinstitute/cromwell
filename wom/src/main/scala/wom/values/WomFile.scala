package wom.values

import java.io.FileNotFoundException
import java.nio.file.NoSuchFileException

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wom.expression.IoFunctionSet
import wom.types._

import scala.concurrent.Future
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
  def withSize(ioFunctionSet: IoFunctionSet): Future[WomFile] = Future.successful(this)

  protected val recoverFileNotFound: PartialFunction[Throwable, this.type] = {
    case _: NoSuchFileException | _: FileNotFoundException => this
    case e if e.getCause != null && recoverFileNotFound.isDefinedAt(e.getCause) => recoverFileNotFound.apply(e.getCause)
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

/**
  * Mix in this trait when creating a WomFile to signal that its actionable path cannot be determined when the object is instantiated.
  * Instead, provide an initialize method that needs to be called when appropriate to perform initialization.
  * The returned value should NOT be a LazyWomFile
  */
trait LazyWomFile { this: WomFile =>
  def initialize(ioFunctionSet: IoFunctionSet): ErrorOr[WomFile]
}

object LazyWomFile {
  implicit class InitializableWomValue(val womValue: WomValue) extends AnyVal {
    def initialize(ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = womValue match {
      case lazyFile: LazyWomFile => lazyFile.initialize(ioFunctionSet)
      case array: WomArray => array.traverse(_.initialize(ioFunctionSet))
      case array: WomObjectLike => array.traverse(_.initialize(ioFunctionSet))
      case map: WomMap => map.traverseValues(_.initialize(ioFunctionSet))
      case optional: WomOptionalValue => optional.traverse(_.initialize(ioFunctionSet))
      case other => other.validNel
    }
  }

  implicit class InitializableWomFile(val womFile: WomFile) extends AnyVal {
    def initializeWomFile(ioFunctionSet: IoFunctionSet): ErrorOr[WomFile] = womFile match {
      case lazyFile: LazyWomFile => lazyFile.initialize(ioFunctionSet)
      case other => other.validNel
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

  override def withSize(ioFunctionSet: IoFunctionSet): Future[WomFile] = {
    ioFunctionSet.size(value)
      .map(s => WomMaybePopulatedFile(valueOption = Option(value), sizeOption = Option(s)))(ioFunctionSet.ec)
      .recover(recoverFileNotFound)(ioFunctionSet.ec)
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
                                   basename: Option[String] = None) extends WomFile {
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
                                 sizeOption: Option[Long] = None,
                                 formatOption: Option[String] = None,
                                 contentsOption: Option[String] = None,
                                 secondaryFiles: Seq[WomFile] = Vector.empty) extends WomFile {
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

  override def withSize(ioFunctionSet: IoFunctionSet): Future[WomMaybePopulatedFile] = {
    implicit val ec = ioFunctionSet.ec
    (sizeOption, contentsOption) match {
      case (Some(_), _) => Future.successful(this)
      case (None, Some(contents)) => Future.successful(this.copy(sizeOption = Option(contents.length.toLong)))
      case _ => ioFunctionSet.size(value)
        .map(s => this.copy(sizeOption = Option(s)))(ec)
        .recover(recoverFileNotFound)(ec)
    }
  }
}

object WomMaybePopulatedFile {
  def apply(value: String): WomMaybePopulatedFile = WomMaybePopulatedFile(valueOption = Option(value))
}
