package perf

import better.files.File
import io.circe
import io.circe.generic.auto._
import io.circe.parser._

//Do not remove this unused import, as removing it fails to find implicit decoder for OffsetDateTime
import io.circe.java8.time.decodeOffsetDateTimeDefault

object CompareMetadata extends App {

  def parseMetadataFromFile(filePath: String): Either[circe.Error, Metadata] = {
//    val metadataFile = File("/Users/sshah/Documents/perf_metadata_compare/small_metadata.json")
    val metadataFile = File(filePath)
    val metadataFileContent = metadataFile.contentAsString

    decode[Metadata](metadataFileContent)
  }

  def computeMetricsFromMetadata(metadata: Metadata): Unit = {
    println(metadata.avgTimeInCallCachingState)
  }

  def printErrorToConsoleAndExit(errorList: circe.Error*): Unit = {
    Console.err.println(s"Something went wrong while parsing metadata.json. Error:")
    errorList.foreach(e => e.printStackTrace(Console.err))
    System.exit(1)
  }


  args.length match {
    //TODO: remove case 1 later
    case 1 => parseMetadataFromFile(args(0))
    case 2 => {
      val metadata1Either = parseMetadataFromFile(args(0))
      val metadata2Either = parseMetadataFromFile(args(1))
      (metadata1Either, metadata2Either) match {
        case (Right(metadata1), Right(metadata2)) => {
          computeMetricsFromMetadata(metadata1)
          computeMetricsFromMetadata(metadata2)
        }
        case (Right(_), Left(e)) => printErrorToConsoleAndExit(e)
        case (Left(e), Right(_)) => printErrorToConsoleAndExit(e)
        case (Left(e1), Left(e2)) => printErrorToConsoleAndExit(e1, e2)
      }
    }
    case _ => {
      Console.err.println("Please pass in only 2 file paths!")
      System.exit(1)
    }
  }
}
