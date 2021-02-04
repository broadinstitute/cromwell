package cromwell.filesystems.gcs.batch

import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import com.google.api.client.googleapis.json.GoogleJsonError
import com.google.api.client.http.HttpHeaders
import com.google.api.services.storage.model.{Objects, RewriteResponse, StorageObject}
import cromwell.filesystems.gcs.{GcsPathBuilder, MockGcsPathBuilder}
import org.scalatest.BeforeAndAfterAll
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.collection.JavaConverters._

class GcsBatchIoCommandSpec extends AnyFlatSpec with Matchers with BeforeAndAfterAll {
  behavior of "GcsBatchIoCommand"

  private lazy val gcsPathBuilder: GcsPathBuilder = MockGcsPathBuilder.instance

  private lazy val gcsPath = gcsPathBuilder.build("gs://hello/world").get

  it should "test GcsBatchDeleteCommand" in {
    val command = PartialGcsBatchCommandBuilder.deleteCommand((gcsPath, false)).get

    type commandType = com.google.api.services.storage.Storage#Objects#Delete
    command.operation should be(a[commandType])

    val operation = command.operation.asInstanceOf[commandType]
    operation.getBucket should be("hello")
    operation.getObject should be("world")

    command.setUserProject should be(false)
    command.withUserProject.setUserProject should be(true)

    command.mapGoogleResponse(null)

    command.onSuccess(null, new HttpHeaders()).toEither.right.get.left.get

    command.onFailure(new GoogleJsonError(), new HttpHeaders()) should be(None)
  }

  it should "test GcsBatchSizeCommand" in {
    val command = PartialGcsBatchCommandBuilder.sizeCommand(gcsPath).get

    type commandType = com.google.api.services.storage.Storage#Objects#Get
    command.operation should be(a[commandType])

    val operation = command.operation.asInstanceOf[commandType]
    operation.getBucket should be("hello")
    operation.getObject should be("world")

    command.setUserProject should be(false)
    command.withUserProject.setUserProject should be(true)

    val response = new StorageObject()
    command.mapGoogleResponse(response) should
      be(Invalid(NonEmptyList.one("'gs://hello/world' in project 'cromwell-test' returned null size")))

    response.setSize(BigInt(139).bigInteger)
    command.mapGoogleResponse(response) should be(Valid(139L))

    command.onSuccess(response, new HttpHeaders()).toEither.right.get.left.get should be(139L)

    command.onFailure(new GoogleJsonError(), new HttpHeaders()) should be(None)
  }

  it should "test GcsBatchCrc32Command" in {
    val command = PartialGcsBatchCommandBuilder.hashCommand(gcsPath).get

    type commandType = com.google.api.services.storage.Storage#Objects#Get
    command.operation should be(a[commandType])

    val operation = command.operation.asInstanceOf[commandType]
    operation.getBucket should be("hello")
    operation.getObject should be("world")

    command.setUserProject should be(false)
    command.withUserProject.setUserProject should be(true)

    val response = new StorageObject()
    response.setCrc32c("aeiouy")

    command.mapGoogleResponse(response) should be(Valid("aeiouy"))

    command.onSuccess(response, new HttpHeaders()).toEither.right.get.left.get should be("aeiouy")

    command.onFailure(new GoogleJsonError(), new HttpHeaders()) should be(None)
  }

  it should "test GcsBatchTouchCommand" in {
    val command = PartialGcsBatchCommandBuilder.touchCommand(gcsPath).get

    type commandType = com.google.api.services.storage.Storage#Objects#Get
    command.operation should be(a[commandType])

    val operation = command.operation.asInstanceOf[commandType]
    operation.getBucket should be("hello")
    operation.getObject should be("world")

    command.setUserProject should be(false)
    command.withUserProject.setUserProject should be(true)

    val response = new StorageObject()
    command.mapGoogleResponse(response)

    command.onSuccess(response, new HttpHeaders()).toEither.right.get.left.get

    command.onFailure(new GoogleJsonError(), new HttpHeaders()) should be(None)
  }

  it should "test GcsBatchExistsCommand" in {
    val command = PartialGcsBatchCommandBuilder.existsCommand(gcsPath).get

    type commandType = com.google.api.services.storage.Storage#Objects#Get
    command.operation should be(a[commandType])

    val operation = command.operation.asInstanceOf[commandType]
    operation.getBucket should be("hello")
    operation.getObject should be("world")

    command.setUserProject should be(false)
    command.withUserProject.setUserProject should be(true)

    val response = new StorageObject()
    command.mapGoogleResponse(response) should be(Valid(true))

    command.onSuccess(response, new HttpHeaders()).toEither.right.get.left.get should be(true)

    val error = new GoogleJsonError()
    command.onFailure(error, new HttpHeaders()) should be(None)

    error.setCode(404)
    command.onFailure(error, new HttpHeaders()) should be(Option(Left(false)))
  }

  it should "test GcsBatchIsDirectoryCommand" in {
    val command = PartialGcsBatchCommandBuilder.isDirectoryCommand(gcsPath).get

    type commandType = com.google.api.services.storage.Storage#Objects#List
    command.operation should be(a[commandType])

    val operation = command.operation.asInstanceOf[commandType]
    operation.getBucket should be("hello")
    operation.getPrefix should be("world/")

    command.setUserProject should be(false)
    command.withUserProject.setUserProject should be(true)

    val response = new Objects()
    command.mapGoogleResponse(response) should be(Valid(false))

    response.setItems(List(new StorageObject()).asJava)
    command.mapGoogleResponse(response) should be(Valid(true))

    command.onSuccess(response, new HttpHeaders()).toEither.right.get.left.get should be(true)

    command.onFailure(new GoogleJsonError(), new HttpHeaders()) should be(None)
  }

  it should "test GcsBatchCopyCommand" in {
    val gcsDestination = gcsPathBuilder.build("gs://howdy/universe").get

    val command = PartialGcsBatchCommandBuilder.copyCommand((gcsPath, gcsDestination)).get

    type commandType = com.google.api.services.storage.Storage#Objects#Rewrite
    command.operation should be(a[commandType])

    val operation = command.operation.asInstanceOf[commandType]
    operation.getSourceBucket should be("hello")
    operation.getSourceObject should be("world")
    operation.getDestinationBucket should be("howdy")
    operation.getDestinationObject should be("universe")

    command.setUserProject should be(false)
    command.withUserProject.setUserProject should be(true)

    val response = new RewriteResponse()
    command.mapGoogleResponse(response) should be(Valid(()))

    response.setDone(true)
    command.onSuccess(response, new HttpHeaders()).toEither.right.get.left.get should be(())

    response.setDone(false)
    response.setRewriteToken("token")
    command.onSuccess(response, new HttpHeaders()).toEither.right.get.right.get.rewriteToken should be(Option("token"))

    command.onFailure(new GoogleJsonError(), new HttpHeaders()) should be(None)
  }

  it should "test ExceptionSpewingGcsBatchIoCommand" in {
    val command = ExceptionSpewingGcsBatchIoCommand

    the[UnsupportedOperationException] thrownBy {
      command.operation
    } should have message "operation is not supported"

    the[UnsupportedOperationException] thrownBy {
      command.withUserProject
    } should have message "withUserProject is not supported"

    the[RuntimeException] thrownBy {
      command.mapGoogleResponse(null)
    } should have message "Ill behaved code that throws in mapGoogleResponse"

    the[RuntimeException] thrownBy {
      command.onSuccess(null, new HttpHeaders())
    } should have message "Ill behaved code that throws in mapGoogleResponse"

    the[RuntimeException] thrownBy {
      command.onFailure(new GoogleJsonError(), new HttpHeaders())
    } should have message "Ill behaved code that throws in onFailure"
  }
}
