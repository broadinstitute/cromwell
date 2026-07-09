AWS/CROMWELL DEPLOYMENT
=======================

Overview
--------

AWS Cloudformation is used to deploy Cromwell-AWS.  As resources are created, costs might be generated along the way.  To make sure all generated data is cleaned, remove the stacks again using the CloudFormation Console.


Deployment
----------

Deployment of the cromwell/AWS environment can be performed using the three cloudformation stacks:

1. VPC : setup of the networks
2. Resources : setup of the compute environment, job queues and storage solutions
3. Cromwell : setup of an EC2 instance and RDS, hosting the cromwell server and submission tools. 

Along the way, all necessary IAM rols are generated. 

*Note:* This tutorial uses the eu-west-2 (London) region. Change where appropriate for you setup

*Note:* Select "Preserve successfully provisioned resources" as failure option to review problems in depth. 

### STEP 1 : SETUP VPC 

The default VPC settings create both private and public subnets. Private subnets are used for the compute environment and database setup. Public subnets are used to host the cromwell server instance. 


Used Template : 

```
https://cromwell-aws-cloudformation-templates.s3.eu-west-1.amazonaws.com/root-templates/aws-vpc.template.yaml
```

Steps:

1. Go to the CloudFormation console in your region : https://eu-west-2.console.aws.amazon.com/cloudformation/home?region=eu-west-2#/

2. Select "Stacks" in the left menu, then "Create Stack" on the top right. Choose "With new resources" :[link](https://eu-west-2.console.aws.amazon.com/cloudformation/home?region=eu-west-2#/stacks/create)

3. Select "Template is ready", and use "Amazon S3 URL" as template source. Provide the yaml file URL provided above.

4. Enter a stack name. 

5. Select availability zones, more is better. This gives Batch more zones to find suitable instances. 

6. Match "Number of availability zones" to your selection. 

7. Other settings can be left as default. 

8. Review, confirm and submit on final page. 



### STEP 2 : SETUP COMPUTE RESOURCES

The default setup configures a spot and on-demand queue, a cromwell-specific bucket and the necessary IAM roles to run BATCH jobs. See below for further configuration to extend this setup.

Used Template : 

```
https://cromwell-aws-cloudformation-templates.s3.eu-west-1.amazonaws.com/root-templates/gwfcore-root.template.yaml
```

Steps:

1. Go to the CloudFormation console in your region : https://eu-west-2.console.aws.amazon.com/cloudformation/home?region=eu-west-2#/

2. Select "Stacks" in the left menu, then "Create Stack" on the top right. Choose "With new resources" :[link](https://eu-west-2.console.aws.amazon.com/cloudformation/home?region=eu-west-2#/stacks/create)

3. Select "Template is ready", and use "Amazon S3 URL" as template source. Provide the yaml file URL provided above.

4. Enter a stack name. 
    * *Note:* This name is referred to as "GWFCORE NameSpace" in STEP 3
    * *Note:* The name gets reflected in queue names etc. Keep it short. 

5. Provide an s3 bucket name to host cromwell temp data, results and runtime requirements. 
    * *Note:* If the bucket exists, it must be located in the same region. 
    * *Note:* If the bucket exists, specify this in the next field

6. Select the VPC ID generated in STEP 1

7. Provide the PRIVATE subnets (all of them) from STEP 1, for the compute environment. Match the number of subnets to your selection

8. Set the Max vCPU count for default (spot) and High Priority (on-demand) compute environment.
    * *Note:*  Check your quota per region [here](https://eu-west-2.console.aws.amazon.com/servicequotas/home/services/ec2/quotas), look for "Standard"

9. Set the maximal spot bidding price.

10. Select the list of instance types to be used. Default "optimal" value is a safe bet. 
    * *Note:*  To edit this, clone the compute environment in the [Batch dashboard](https://eu-west-2.console.aws.amazon.com/batch/home?region=eu-west-2#compute-environments), and play around with the options. 

11. *Optional:* Create an EFS filesystem:
    * *Note:* EFS is a distributed filesystem. Keeping large intermediate results on EFS can improve perforance by reducing S3/EC2 transfer times
    * *Note:* By default, EFS performance is limited. Change the setup of the volume using the [console](https://eu-west-2.console.aws.amazon.com/efs/home?region=eu-west-2#/file-systems) to "Enhanced/Elastic" Performance. Consider the costs this implies! 
    * *Note:* It's recommended to create a new EFS, to make sure that mounting and network settings are correctly setup. 


12. *Optional:* Create an FSx filesystem:
    * *Note:* see documentation in [README](README.md)

13. Other settings can be left as default. 

14. Review, confirm and submit on final page. 


### STEP 3 : SETUP CROMWELL SERVER

The default setup deploys a small t3.medium EC2 instance with 25Gb of storage and a RDS aurora-mysql database, to host cromwell, and some interaction tools. Although enough for testing, you'll probably have to scale up the instance for production runs. Alternatively, you can deploy cromwell on local infrastructure, preventing the RDS/EC2 costs.  See below for details.


Used Template : 

```
https://cromwell-aws-cloudformation-templates.s3.eu-west-1.amazonaws.com/root-templates/cromwell-resources.template.yaml
```

Steps:

1. Go to the CloudFormation console in your region : https://eu-west-2.console.aws.amazon.com/cloudformation/home?region=eu-west-2#/

2. Select "Stacks" in the left menu, then "Create Stack" on the top right. Choose "With new resources" :[link](https://eu-west-2.console.aws.amazon.com/cloudformation/home?region=eu-west-2#/stacks/create)

3. Select "Template is ready", and use "Amazon S3 URL" as template source. Provide the yaml file URL provided above.

4. Provide a stack name.

5. Provide a name space.
    * *Note:* This name is reflected in the key name. Keep it short.

6. Provide the GWFCORE Namespace from STEP 2

7. Specify the VPC ID, from STEP 1

8. Select a PUBLIC subnet for the cromwell server

9. Select 2+ PRIVATE subnets for the RDS database

10. Set the instance root volume size to 25+ Gb

11. Set a password for the cromwell database

12. Provide Filesystem details for EFS if specified in STEP 2:
    * *Note:* Get the values from [the dashboard](https://eu-west-2.console.aws.amazon.com/efs/home?region=eu-west-2#/file-systems)
    * *Note:* Accesspoint is listed under filesystem details, tab "Access Points"
    * *Note:* Security Group is listed under filesystem details, tab "Network"

13. Other settings can be left as default. 

14. Review, confirm and submit on final page. 


When the stack is ready, use the following commands to retrieve your SSH key. The exact name of the key can be retrieved from the EC2 instance "Connect" page. 

```
KEY_NAME=key-<stack_namespace>
REGION=<region>

# GET KEY ID :
KEY_ID=$(aws ec2 describe-key-pairs --filters Name=key-name,Values=${KEY_NAME} --query KeyPairs[*].KeyPairId --output text --region ${REGION})

# GET KEY CONTENTS
mkdir -p ~/.ssh/.keys
aws ssm get-parameter --region ${REGION} --name /ec2/keypair/${KEY_ID} --with-decryption --query Parameter.Value --output text > ~/.ssh/.keys/${KEY_NAME}.pem
chmod 600 ~/.ssh/.keys/${KEY_NAME}.pem

```



Post Install Optimizations
--------------------------

### EFS 

EFS data is not cleaned up automatically. Consider adding a "cleanup" step in your production WFs, to keep the storage footprint/costs low.

Edit settings [here](https://eu-west-2.console.aws.amazon.com/efs/home?region=eu-west-2#/file-systems). Select the FS and "edit"

* Performance : "Enhanced":"Elastic" Throughput mode is recommended for large workflows. consider extra costs
* Lifecycle management : Consider moving data to infrequent access, to reduce costs of "forgotten" data. 

### S3 

The cromwell temporary data and results are located in the bucket specified in STEP 2. To reduce your storage footprint/costs, consider setting up lifecycle management. 

* *cromwell-execution/* : Contains analysis results and logs for each executed workflow. 
* *scripts/* : contains the runtime scripts for each executed task 

To set up a lifecycle:

* go to the [s3 console](https://s3.console.aws.amazon.com/s3/buckets?region=eu-west-2), and select the bucket from STEP 1
* Select the "Management" tab and "Create Lifecycle rule"
* analysis results:
    * Set name : e.g. "Prune cromwell-execution"
    * Set prefix "cromwell-execution"
    * Choose actions : "Expire current versions of objects" and "Delete expired object delete markers or incomplete multipart uploads"
    * Set time limits for  after which objects are deleted (expired).
    * Enable "delete incomplete multipart upload" checkbox and set time. 
* scripts:
    * repeat steps above, with prefix "scripts"


### LOCAL CROMWELL INSTANCE

Running cromwell locally has the benefit of not running an EC2 instance 24/7.  However, consider the following points: 

* Run cromwell as a user with a default AWS IAM profile set to the same region, and with sufficient permissions to interact with ECR, Batch and S3. 

* Setup a local database, or setup ingress rules on the RDS security group. 

* Consider load balancing by running multiple instances : see [here](https://cromwell.readthedocs.io/en/stable/Scaling/)

* for EFS, there are specific concerns:
    * expose EFS to your local network : 
        * Run a cheap (e.g t4g.nano) instance mounting the EFS share, assign a PEM key to access it
        * Setup sshfs to mount the EFS share from that instance to your local machine at "/mnt/efs"
    * REDUCE network traffic by enabling the "check-sibling-md5" and "efsMakeMD5" settings


### COMPUTE ENVIRONMENT

* You might play around with the instance types to get more options:
    * clone the compute environment
    * double check all network settings ! 
    * replace the "optimal" type by "*.family" types. There is maximal number of entries, so using whole families allows more types

* Extra Queues : You might consider a dedicated queue for high disk/network jobs, by altering the launch template:
    * Go to Batch console, select the spot compute environment and open the JSON tab. Look for "launchTemplateName" (eg lt-06fa9fee031254098)
    * Go to EC2 console, select "Launch Templates" in the menu on the left, and search for your launch template
    * On details, select "Modify Template (create new version)" in the top right under "Actions"
    * On the resulting page, add a description (eg "HighNetwork"), and then open advanced settings at the bottom of the page
    * At the bottom, look in "User Data". About halfway, set EBS_* settings. EBS_IOPS can go as high as 16,000 and EBS_THROUGHPUT can go as high as 1,000Mbs/s
    * Now clone the compute environment (double check all network settings and roles)
        * specify "exact" launch template version to the version you created. 
        * set instance type to network/disk optimized machines (eg m5zn.6xlarge)
        * have a blazing fast I/O machine (tests reached constant simultaneous 1Gb/s upload and 1gb/s download, while transferring novaseq data from basespace to AWS/S3)


 


