package centaur.test

import java.io.File

sealed trait FilesChecker {
  def countObjectsAtPath: String => Int
}

object LocalFilesChecker extends FilesChecker {
  def countObjectsAtPath: String => Int = { s =>
    val d = new File(s)
    if (d.exists && d.isDirectory)
      d.listFiles.length
    else if (d.exists && d.isFile)
      1
    else
      0
  }
}

object PipelinesFilesChecker extends FilesChecker {

  import ObjectCounterInstances.gcsObjectCounter
  import ObjectCounterSyntax._

  private lazy val storage = Operations.storage

  private val gsPrefixRegex = "^gs:\\/\\/.*"

  def countObjectsAtPath: String => Int = ObjectCounterSyntax(storage).countObjects(gsPrefixRegex)
}

object AWSFilesChecker extends FilesChecker {

  import ObjectCounterInstances.awsS3ObjectCounter
  import ObjectCounterSyntax._

  private lazy val s3Client = Operations.buildAmazonS3Client

  private val s3PrefixRegex = "^s3:\\/\\/.*"

  override def countObjectsAtPath: String => Int = s3Client.countObjects(s3PrefixRegex)
}
