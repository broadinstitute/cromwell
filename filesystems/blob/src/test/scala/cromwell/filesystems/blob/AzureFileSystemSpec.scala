package com.azure.storage.blob.nio;

import cromwell.filesystems.blob._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.nio.file.spi.FileSystemProvider
import java.time.Instant
import scala.compat.java8.OptionConverters._
import scala.jdk.CollectionConverters._

class AzureFileSystemSpec extends AnyFlatSpec with Matchers {
  val now = Instant.now()
  val container = BlobContainerName("testConainer")
  val exampleSas = BlobPathBuilderFactorySpec.buildExampleSasToken(now)
  val exampleConfig = BlobFileSystemManager.buildConfigMap(exampleSas, container)
  val exampleEndpoint = BlobPathBuilderSpec.buildEndpoint("testStorageAccount")
  val fileSystemProvider = FileSystemProvider.installedProviders().asScala.find(p => p.getScheme() == "azb")
  it should "parse an expiration from a sas token" in {
    fileSystemProvider.nonEmpty shouldBe(true)
    FileSystemProvider.installedProviders().asScala.map(_.getScheme) shouldBe List("azb")
    val fs = fileSystemProvider.map(p => new AzureFileSystem(p.asInstanceOf[AzureFileSystemProvider], exampleEndpoint.value, exampleConfig.asJava))
    fs.nonEmpty shouldBe(true)
    fs.flatMap(_.getExpiry().asScala) shouldBe(Some(now))
    fs.map(_.getFileStore().name()) shouldBe(Some(container.value))
  }
}
