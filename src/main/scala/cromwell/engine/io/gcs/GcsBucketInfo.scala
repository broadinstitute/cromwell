package cromwell.engine.io.gcs

import com.google.api.client.util.DateTime
import com.google.api.services.storage.model.Bucket.Owner

case class GcsBucketInfo(bucketName: String, location: String, timeCreated: DateTime, owner: Owner)