package wom.values

import java.io.FileNotFoundException
import java.nio.file.NoSuchFileException

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
    * Returns the WomPrimitiveFile instances recursively referenced by this instance.
    *
    * WomMaybeListedDirectory instances return either just the directory as an WomUnlistedDirectory, or if there is a
    * listing then returns the recursive listing. WomMaybePopulatedFile instances return that instance as a
    * WomSingleFile plus any primitives recursively discovered in secondaryFiles.
    *
    * WomPrimitiveFile instances return just the instance.
    */
  def flattenFiles: Seq[WomPrimitiveFile] =
    this match {
      case womPrimitiveFile: WomPrimitiveFile => List(womPrimitiveFile)
    }

  def isFileNotFound(t: Throwable): Boolean = t match {
    case _: NoSuchFileException | _: FileNotFoundException => true
    case other if other.getCause != null => isFileNotFound(t.getCause)
    case _ => false
  }
}

object WomFile {
  def apply(fileType: WomFileType, value: String) =
    fileType match {
      case WomUnlistedDirectoryType => WomUnlistedDirectory(value)
      case WomSingleFileType => WomSingleFile(value)
      case WomGlobFileType => WomGlobFile(value)
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
    case r: WomString => Success(this.copy(value = value + r.value.trim))
    case r: WomOptionalValue => evaluateIfDefined("+", r, add)
    case _ => invalid(s"$value + $rhs")
  }

  override def equals(rhs: WomValue): Try[WomBoolean] = rhs match {
    case r: WomUnlistedDirectory => Success(WomBoolean(this.equals(r)))
    case r: WomString => Success(WomBoolean(value.equals(r.value.trim)))
    case r: WomOptionalValue => evaluateIfDefined("==", r, equals)
    case _ => invalid(s"$value == $rhs")
  }

  override def mapFile(f: String => String): WomUnlistedDirectory =
    this.copy(value = f(value))

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
    case r: WomString => Success(this.copy(value = value + r.value.trim))
    case r: WomOptionalValue => evaluateIfDefined("+", r, add)
    case _ => invalid(s"$value + $rhs")
  }

  override def equals(rhs: WomValue): Try[WomBoolean] = rhs match {
    case r: WomSingleFile => Success(WomBoolean(this.equals(r)))
    case r: WomString => Success(WomBoolean(value.equals(r.value.trim)))
    case r: WomOptionalValue => evaluateIfDefined("==", r, equals)
    case _ => invalid(s"$value == $rhs")
  }

  override def mapFile(f: String => String): WomSingleFile =
    this.copy(value = f(value))

  override def mapWomFile(f: WomFile => String) =
    this.copy(value = f(this))
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
    case r: WomString => Success(this.copy(value + r.value.trim))
    case r: WomOptionalValue => evaluateIfDefined("+", r, add)
    case _ => invalid(s"$value + $rhs")
  }

  override def equals(rhs: WomValue): Try[WomBoolean] = rhs match {
    case r: WomGlobFile => Success(WomBoolean(value.equals(r.value) && womType.equals(r.womType)))
    case r: WomString => Success(WomBoolean(value.toString.equals(r.value.toString.trim)))
    case r: WomOptionalValue => evaluateIfDefined("==", r, equals)
    case _ => invalid(s"$value == $rhs")
  }

  override def mapFile(f: String => String): WomGlobFile = this.copy(value = f(value))

  override def mapWomFile(f: WomFile => String) = this.copy(value = f(this))
}
