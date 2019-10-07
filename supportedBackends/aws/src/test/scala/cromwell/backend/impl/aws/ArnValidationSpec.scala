package cromwell.backend.impl.aws

import cats.data.NonEmptyList
import cats.data.Validated.{Invalid, Valid}
import org.scalatest.{Matchers, WordSpecLike}
import wom.values.WomString

class ArnValidationSpec extends WordSpecLike with Matchers{
  private val arnKey = "arn"
  private val arnValidator = ArnValidation(arnKey)

  "ArnValidation" should {

    "validate a valid arn entry" in {
      val validArnsAsStrings = List(
        "arn:aws:batch:us-east-1:111122223333:job-queue/HighPriority",
        "arn:aws:batch:us-west-2:123456789012:job-queue/default-a4e50e00-b850-11e9",
        "arn:aws:batch:us-west-2:123456789012:job-queue:default-a4e50e00-b850-11e9",
        "arn:aws:batch:us-west-2:123456789012:job-queue:default-a4e50e00-b850-11e9:1",
        "arn:aws:batch:ap-northeast-2:123456789012:job-queue:default-a4e50e00-b850-11e9",
        "arn:aws:batch:us-gov-west-1:123456789012:job-queue/default-a4e50e00-b850-11e9",
        "arn:aws:batch:us-gov-west-1:123456789012:job-queue:default-a4e50e00-b850-11e9",
        "arn:aws-cn:batch:us-west-2:123456789012:job-queue/default-a4e50e00-b850-11e9",
        "arn:aws-cn:batch:us-west-2:123456789012:job-queue:default-a4e50e00-b850-11e9",
        "arn:aws-cn:batch:us-west-2:123456789012:default-a4e50e00-b850-11e9",
        "arn:aws-cn:batch:us-gov-west-1:123456789012:job-queue/default-a4e50e00-b850-11e9",
        "arn:aws-cn:batch:us-gov-west-1:123456789012:job-queue:default-a4e50e00-b850-11e9",
        "arn:aws-us-gov:batch:us-west-2:123456789012:job-queue/default-a4e50e00-b850-11e9",
        "arn:aws:batch:us-east-1:123456789012:compute-environment/my-environment",
        "arn:aws:batch:us-east-1:123456789012:job-definition/my-job-definition:1",
        "arn:aws:batch:us-east-1:123456789012:job-queue/my-queue",
        "arn:aws:s3:::my_corporate_bucket/exampleobject.png",
        "arn:aws:elasticbeanstalk:us-east-1:123456789012:environment/My App/MyEnvironment",
        "arn:aws:iam::123456789012:user/David",
        "arn:aws:rds:eu-west-1:123456789012:db:mysql-db",
        "arn:aws:a4b:us-east-1:123456789012:room/7315ffdf0eeb874dc4ab8a546e8b70ec/5f90e5d608b6baa9c88db56654aef158",
        "arn:aws:ecs:us-east-1:123456789012:container-instance/my-cluster/403125b0-555c-4473-86b5-65982db28a6d",
        "arn:aws:resource-groups:us-west-2:123456789012:group/MyExampleGroup/Myfile.png",
        "arn:aws:s3:::my_corporate_bucket/*",
        "arn:aws:s3:::my_corporate_bucket/Development/*",
        "arn:aws:autoscaling:us-east-1:123456789012:scalingPolicy:c7a27f55-d35e-4153-b044-8ca9155fc467:autoScalingGroupName/my-test-asg1:policyName/my-scaleout-policy",
        "arn:aws:waf-regional:us-east-1:123456789012:rule/41b5b052-1e4a-426b-8149-3595be6342c2",
        "arn:aws:lambda:us-west-2:123456789012:function:helloworld:$LATEST",
        "arn:aws:swf:us-east-1:123456789012:/domain/department1",
        "arn:aws:swf:*:123456789012:/domain/*",
        "arn:aws:cognito-sync:us-east-1:123456789012:identitypool/us-east-1:1a1a1a1a-ffff-1111-9999-12345678"
      )
      validArnsAsStrings foreach { arn =>
        val keyToValue = Map(arnKey -> WomString(arn))
        arnValidator.validate(keyToValue) shouldBe Valid(arn)
      }
    }

    "fail to validate an invalid arn entry" in {
      val invalidArnsAsStrings = List(
        "arn:aws:iam::123456789012:u*",
        "arn:aws:s3::my_corporate_bucket",
        "arn:AWS:batch:us-west-2:123456789012:job-queue/default-a4e50e00-b850-11e9",
        "arn:aws-CN:batch:us-west-2:123456789012:job-queue/default-a4e50e00-b850-11e9",
        "arn:aws:waf_regional:us-east-1:123456789012:rule/41b5b052-1e4a-426b-8149-3595be6342c2",
        "arn:aws:s3:::my_corporate_bucket/*text"
      )
      invalidArnsAsStrings foreach { arn =>
        val keyToValue = Map(arnKey -> WomString(arn))
        arnValidator.validate(keyToValue) shouldBe Invalid(NonEmptyList("ARN has invalid format", Nil))
      }
    }
  }
}
