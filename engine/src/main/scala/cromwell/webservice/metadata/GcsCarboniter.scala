package cromwell.webservice.metadata

import akka.actor.ActorRef
import cats.data.{Chain, IndexedReaderWriterStateT}
import cats.effect.IO
import com.google.cloud.storage.{Blob, BlobId, BlobInfo, Storage}


/**
  * State machine for carboniting a workflow from
  */
object GcsCarboniter {

  //This is a process, defined by types of state.
  case class Uuid(value: String)
  case class MetadataRetrieved(uuid: String, json: String)
  case class StoredInGcs(uuid: String, blob: Blob)

  case class Config(storage: Storage, bucket: String, metadataActor: ActorRef)

  def metadata: IndexedReaderWriterStateT[IO, Config, Chain[String], Uuid, MetadataRetrieved, Unit] =
    IndexedReaderWriterStateT[IO, Config, Chain[String], Uuid, MetadataRetrieved, Unit] ( (config, uuid) => {
      val _uuid: String = uuid.value
      //config.metadataActor ? GetMeMetadata(_uuid)
      val metadata: String = null

      IO((Chain.empty[String], MetadataRetrieved(_uuid, metadata), ()))
    })

  def carbonite: IndexedReaderWriterStateT[IO, Config, Chain[String], MetadataRetrieved, StoredInGcs, Unit] =
    IndexedReaderWriterStateT[IO, Config, Chain[String], MetadataRetrieved, StoredInGcs, Unit] ((config, metadataRetrieved) => {
      val blobId = BlobId.of(config.bucket, s"${metadataRetrieved.uuid}.json")
      val blobInfo = BlobInfo.newBuilder(blobId).setContentType("text/plain").build
      val blob: IO[Blob] = IO { config.storage.create(blobInfo, metadataRetrieved.json.getBytes())}
      blob.map(blob => (Chain.empty[String], StoredInGcs(metadataRetrieved.uuid, blob), ()))
    })
}
