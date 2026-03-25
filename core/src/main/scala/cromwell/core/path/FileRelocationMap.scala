package cromwell.core.path

/**
 * Defines source/destination relationships for files being copied with the final outputs location feature
 * @param map all files to be copied for this workflow
 */
case class FileRelocationMap(map: Map[Path, Path]) {

  /**
   * We do NOT want to recompute this every time (CTM-362)
   */
  lazy val stringified: Map[String, String] = map map { case (src: Path, dst: Path) =>
    src.pathAsString -> dst.pathAsString
  }
}

object FileRelocationMap {
  lazy val empty: FileRelocationMap = FileRelocationMap(Map.empty)
}
