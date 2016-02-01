package cromwell.caching

import java.nio.file.Path

/**
  * Provides basic definition to work with files and container images hashes.
  */
trait Cacheable {
  /**
    * Get file content hash.
    * @param files List of files to compute.
    * @return List of hashes related to files.
    */
  def computeInputFileHash(files: List[Path]): Map[Path, String]

  /**
    * Get container image hash.
    * @param imageName Image name.
    * @return A hash.
    */
  def computeContainerImageHash(imageName: String): Map[String, String]
}
