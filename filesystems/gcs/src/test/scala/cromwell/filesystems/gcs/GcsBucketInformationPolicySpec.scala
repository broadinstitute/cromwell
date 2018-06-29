package cromwell.filesystems.gcs

import com.google.cloud.storage.Storage.{BucketField, BucketGetOption}
import com.google.cloud.storage.{Bucket, Storage}
import com.google.common.cache.CacheBuilder
import cromwell.filesystems.gcs.bucket.GcsBucketInformationPolicies.{CachedPolicy, DisabledPolicy, OnDemandPolicy}
import cromwell.filesystems.gcs.bucket.{GcsBucketInformation, RequesterPaysCached, RequesterPaysValue}
import cromwell.filesystems.gcs.cache.GcsBucketCache.BucketGetter
import org.scalamock.scalatest.MockFactory
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.Mockito

class GcsBucketInformationPolicySpec extends FlatSpec with Matchers with Mockito with MockFactory {

  "DisabledPolicy" should "always return Disabled requester pays value" in {
    DisabledPolicy.toRequestHandler(any[Storage], anyString).requesterPays(anyString) shouldBe RequesterPaysValue.Disabled
  }

  "OnDemandPolicy" should "always return Unknown requester pays value" in {
    OnDemandPolicy.toRequestHandler(any[Storage], anyString).requesterPays(anyString) shouldBe RequesterPaysValue.Unknown
  }

  // We can't easily mock out the BucketGetter function from here, instead make sure we get a RequesterPaysCached and separately unit test it
  "CachedPolicy" should "return a RequesterPaysCached object" in {
    CachedPolicy(CacheBuilder.newBuilder().build()).toRequestHandler(any[Storage], anyString) shouldBe a[RequesterPaysCached]
  }

  "RequesterPaysCached" should "cache requester pays values" in {
    val guavaCache = CacheBuilder.newBuilder().build[String, GcsBucketInformation]()

    //    def buildBucketWithRp = new Bucket.Builder().setRequesterPays(true).build()
    //    def buildBucketWithoutRp = new Bucket.Builder().setRequesterPays(false).build()
    //    def throwNotEnoughPermission = throw new StorageException(DoesNotHaveServiceUsePermissionErrorCode, DoesNotHaveServiceUsePermissionErrorMessage)
    //    def throwBucketHasRequesterPays = throw new StorageException(DoesNotHaveServiceUsePermissionErrorCode, DoesNotHaveServiceUsePermissionErrorMessage)
    //    
    val mockedBucketGetter = mockFunction[(String, List[BucketGetOption]), Bucket]
    val projectAndBillingOption = List(BucketGetOption.userProject("project-id"), BucketGetOption.fields(BucketField.BILLING))


    val cache = new RequesterPaysCached(any[Storage], "project-id", guavaCache) {
      override lazy val bucketGetter: BucketGetter = Function.untupled(mockedBucketGetter)
    }
    
    mockedBucketGetter.expects(("rp-bucket", projectAndBillingOption)).returns(new Bucket.Builder().setRequesterPays(true).build())
    //    val cache = new RequesterPaysCached(any[Storage], "project-id", guavaCache) {
    //      override lazy val bucketGetter: BucketGetter = {
    //        // Bucket has RP, and we have permission to get metadata
    //        case ("rp-bucket", options) if
    //        options.contains(BucketGetOption.userProject("project-id")) &&
    //          options.contains(BucketGetOption.fields(BucketField.BILLING)) => buildBucketWithRp
    //
    //        // Bucket does not have RP, and we have permission to get metadata  
    //        case ("non-rp-bucket", options) if
    //        options.contains(BucketGetOption.userProject("project-id")) &&
    //          options.contains(BucketGetOption.fields(BucketField.BILLING)) => buildBucketWithoutRp
    //
    //        // We don't know if bucket has RP because we don't have permission  
    //        case ("not-enough-permission-rp" | "not-enough-permission-no-rp", options) if
    //        options.contains(BucketGetOption.userProject("project-id")) &&
    //          options.contains(BucketGetOption.fields(BucketField.BILLING)) =>
    //          throwNotEnoughPermission
    //
    //        // Second attempt to find out if bucket has RP (it does), without project this time  
    //        case ("not-enough-permission-rp", options) if
    //        !options.contains(BucketGetOption.userProject("project-id")) &&
    //          options.contains(BucketGetOption.fields(BucketField.ID)) =>
    //          throwBucketHasRequesterPays
    //
    //        // Second attempt to find out if bucket has RP (it does not), without project this time  
    //        case ("not-enough-permission-no-rp", options) if
    //        !options.contains(BucketGetOption.userProject("project-id")) &&
    //          options.contains(BucketGetOption.fields(BucketField.ID)) =>
    //          buildBucketWithoutRp
    //      }
  }

}
