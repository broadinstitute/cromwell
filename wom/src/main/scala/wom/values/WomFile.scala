package wom.values

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
    * will traverse any files definied within file, also converting the values within those files too.
    *
    * @param f The function to update the location.
    * @return A new WomFile with the updated location.
    * @see [[wom.values.WomValue.collectAsSeq]]
    * @see [[wom.WomFileMapper.mapWomFiles]]
    */
  def mapFile(f: String => String): WomFile

  /**
    * Returns the WomFile recursively referenced by this instance.
    *
    * WomSingleDirectory return either just the directory, or if there is Some listing then just the recursive listing
    * entries. WomSingleFile return that instance plus any instances recursively discovered in secondaryFiles.
    * WomGlobFile return just the instance.
    */
  def flattenFiles: Seq[WomFile] = {
    this match {
      case womMaybeListedDirectory: WomMaybeListedDirectory =>
        womMaybeListedDirectory.listingOption.getOrElse(Nil).toList match {
          case Nil => List(this)
          case list => list.flatMap(_.flattenFiles)
        }
      case womMaybePopulatedFile: WomMaybePopulatedFile =>
        womMaybePopulatedFile.secondaryFiles.foldLeft(List(this)) {
          (womFiles, womMaybeListedDirectoryOrFile) =>
            womFiles ++ womMaybeListedDirectoryOrFile.flattenFiles
        }
      case womPrimitiveFile: WomPrimitiveFile => List(womPrimitiveFile)
    }
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
}


/**
  * A directory possibly with a listing of other files/directories.
  *
  * @param valueOption   The location of the directory, possibly in the cloud.
  * @param listingOption An optional listing of files/directories, either supplied by a user or generated by the engine.
  */
final case class WomMaybeListedDirectory(valueOption: Option[String] = None,
                                         listingOption: Option[Seq[WomFile]] = None) extends WomFile {
  override def value: String = {
    valueOption.getOrElse(throw new UnsupportedOperationException(s"value is not available: $this"))
  }

  override val womFileType: WomFileType = WomMaybeListedDirectoryType

  // TODO: WOM: WOMFILE: This isn't even close to a WDL representation (and maybe belongs in WDL?) of this class, but w/o it other areas of the code crash
  override def toWomString = s""""$value""""

  override def mapFile(f: String => String): WomMaybeListedDirectory = {
    this.copy(valueOption = valueOption.map(f), listingOption.map(_.map(_.mapFile(f))))
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
final case class WomMaybePopulatedFile(valueOption: Option[String] = None,
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
}

object WomMaybePopulatedFile {
  def apply(value: String): WomMaybePopulatedFile = WomMaybePopulatedFile(valueOption = Option(value))
}
