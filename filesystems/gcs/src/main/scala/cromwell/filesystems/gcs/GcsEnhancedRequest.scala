package cromwell.filesystems.gcs

import akka.http.scaladsl.model.StatusCodes
import cats.effect.IO
import com.google.api.client.googleapis.json.GoogleJsonResponseException
import com.google.cloud.storage.StorageException
import cromwell.filesystems.gcs.RequesterPaysErrors.isProjectNotProvidedError

import java.io.FileNotFoundException

object GcsEnhancedRequest {

  // If the request fails because no project was passed, recover the request, this time setting the project
  def recoverFromProjectNotProvided[A](path: GcsPath, f: Boolean => A) = {
    IO(f(false)).handleErrorWith({
        // Only retry with the the project if the error is right
        case error: StorageException if isProjectNotProvidedError(error) =>
          IO(f(true))
        // Use NoSuchFileException for better error reporting
        case e: StorageException if e.getCode == StatusCodes.NotFound.intValue =>
          IO.raiseError(new FileNotFoundException(s"File not found: ${path.pathAsString}"))
        case e: GoogleJsonResponseException if isProjectNotProvidedError(e) =>
          IO(f(true))
        case e: GoogleJsonResponseException if e.getStatusCode == StatusCodes.NotFound.intValue =>
          IO.raiseError(new FileNotFoundException(s"File not found: ${path.pathAsString}"))
        case e =>
          IO.raiseError(e)
      })
  }
}
