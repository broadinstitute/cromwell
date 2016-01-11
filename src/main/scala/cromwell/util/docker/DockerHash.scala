package cromwell.util.docker

import javax.xml.bind.DatatypeConverter

import cromwell.engine.ErrorOr
import cromwell.engine.Hashing._
import cromwell.util.TryUtil

import scala.util.{Failure, Success, Try}
import scalaz.Scalaz._
import scalaz.{Failure => FailureZ, Success => SuccessZ, _}

case class DockerHash(hashType: String, hashString: String) {
  val digest = s"$hashType:$hashString"
}

object DockerHash {
  /**
    * Creates a unique hash from a sequence of hashes.
    * If there is a single hash, returns just the single hash string.
    * If there is more than one hash, concatenates the hash strings together, and then md5s the resulting string.
    * Returns a failure if there are zero hashes, or if there are more than one type of hash type.
    *
    * @param hashCollectionType A description of the collection.
    * @param dockerHashes The sequence of hashes to rehash.
    * @return A single docker hash.
    */
  def fromSeq(hashCollectionType: String, dockerHashes: Seq[DockerHash]): Try[DockerHash] = {
    dockerHashes.size match {
      case 0 =>
        // Need at least one docker hash to hash.
        Failure(new IllegalArgumentException("docker hashes is empty"))

      case 1 =>
        // Just use the original hash, but precede it with the collection type.
        val dockerHash = dockerHashes.head
        Success(dockerHash.copy(s"$hashCollectionType-${dockerHash.hashType}"))

      case _ =>
        // Concatenate the hash strings together, then md5 the result to create another hash.
        val hashTypes = dockerHashes.map(_.hashType).distinct
        hashTypes.size match {
          case 1 =>
            val hashType = s"$hashCollectionType-${hashTypes.head}-md5"
            val hashString = dockerHashes.map(_.hashString).mkString("").md5Sum
            Success(DockerHash(hashType, hashString))
          case _ => Failure(new IllegalArgumentException(s"found more than one docker hash type: $hashTypes"))
        }
    }
  }

  /**
    * Creates a unique hash from a sequence of hashes.
    * If there is a single hash, returns just the single hash string.
    * If there is more than one hash, concatenates the hash strings together, and then md5s the resulting string.
    * Returns a failure if there are zero hashes, or if there are more than one type of hash type.
    *
    * @param hashCollectionType A description of the collection.
    * @param dockerHashes The sequence of hashes to rehash.
    * @return A single docker hash.
    */
  def fromTries(hashCollectionType: String, dockerHashes: Seq[Try[DockerHash]]): Try[DockerHash] = {
    TryUtil.sequence(dockerHashes).flatMap(fromSeq(hashCollectionType, _))
  }

  /**
    * Parses a string like "type:12345678" into a DockerHash(type, hash)
    */
  def fromDigest(digest: String): Try[DockerHash] = {
    digest.indexOf(':') match {
      case -1 => fromHash("unknown", digest)
      case index => fromHash(digest.substring(0, index), digest.substring(index + 1))
    }
  }

  /** Creates a docker hash, after validating the hash string. */
  def fromHash(hashType: String, hashString: String): Try[DockerHash] = {
    validateHashString(hashString) match {
      case SuccessZ(_) => Success(DockerHash(hashType, hashString))
      case FailureZ(e) =>
        val errorMessages = e.toList.mkString(", ")
        Failure(new IllegalArgumentException(s"hashString '$hashString' is not valid: $errorMessages"))
    }
  }

  /** Validates the hash string. */
  private def validateHashString(hashString: String): ErrorOr[String] = {
    val validation = validateHashStringHex(hashString) +++ validateHashStringLength(hashString)
    // Turn the concatenated results back into just the hashString.
    validation map { _ => hashString }
  }

  /** Return the hash if it's valid hex, or the exception message. */
  private def validateHashStringHex(hashString: String): ErrorOr[String] = {
    // We only want to know that we _could_ parse the hash.
    Validation
      .fromTryCatchNonFatal(DatatypeConverter.parseHexBinary(hashString))
      .map(_ => hashString)
      .leftMap(_.getMessage)
      .toValidationNel
  }

  /** Return the hash if it has a valid length, or an error message. */
  private def validateHashStringLength(hashString: String): ErrorOr[String] = {
    val length = hashString.length
    if (length == 8 || length == 32 || length == 64) {
      hashString.successNel
    } else {
      s"unexpected hash length: $length".failureNel
    }
  }
}
