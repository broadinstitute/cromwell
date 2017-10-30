**Spark Backend**

This backend adds support for execution of spark jobs in a workflow.

It supports the following Spark deploy modes:

*  Client deploy mode using the spark standalone cluster manager
*  Cluster deploy mode using the spark standalone cluster manager
*  Client deploy mode using Yarn resource manager
*  Cluster deploy mode using Yarn resource manager

**Configuring Spark Project**

Cromwell's default configuration file is located at `core/src/main/resources/reference.conf`

To customize configuration it is recommended that one copies relevant stanzas from `core/src/main/resources/reference.conf` into a new file, modify it as appropriate, then pass it to Cromwell via:

`java -Dconfig.file=/path/to/yourOverrides.conf cromwell.jar`

Spark configuration stanza is as follows: 

```conf
Spark {
       actor-factory = "cromwell.backend.impl.spark.SparkBackendFactory"
       config {
         # Root directory where Cromwell writes job results.  This directory must be
         # visible and writeable by the Cromwell process as well as the jobs that Cromwell
         # launches.
         root: "cromwell-executions"

         filesystems {
           local {
             localization: [
               "hard-link", "soft-link", "copy"
             ]
           }
         }
		master: "local"
		deployMode: "client"
        }

      }
```
and add backend provider as Spark. 

```
backend {
  default = "Spark"
  providers {
  ....
```

**Configuring Spark Master and Deploy Mode**

Default configuration is as follows:

```conf
Spark {
		......
		master: "local"
		deployMode: "client"

      }
```

However to use Spark in standalone cluster mode change `master: spark://hostname:6066` and `deployMode: cluster` similarly, for yarn change `master: yarn` and `deployMode: cluster` or `deployMode: client` to run in cluster or client mode respectively. 

**Spark runtime attributes**

Supported runtime attributes for a Spark Job is as follows:
	
* executorCores (default value is 1)
* executorMemory (default value is "1GB", `Unit in MB or GB or TB.. ` )
* appMainClass ( Spark app/job entry point)
* numberOfExecutors ( Specific to cluster deploy mode)
* additionalArgs ( i.e to add additional configuration or parameters to spark-submit)

Sample usage:

```wdl
task sparkjob_with_yarn_cluster {
        .....
        
        runtime {
                appMainClass: "${entry_point}"
                executorMemory: "4GB"
                executorCores: "2"
                additionalArgs: "--conf '-Dsamjdk.compression_level=1 -Dsnappy.disable=true' ...."
        }
        
        .....
	}
```

**Spark Environment**

The Spark backend assumes Spark is already installed, and it constructs the spark submit command with the `SPARK_HOME` environment variable if set. Otherwise backend creates command `spark-submit` without a fully qualified path to `spark-submit`.Also, it is important to set environment variable `HOSTNAME` to master machine ip or hostname, that is accessible by spark backend. That can be done by setting either in `~/.bashrc or profile like "export HOSTNAME=<machine ip>" `

Supported File Systems as follows: 

* Local File System
* Network File System
* Distributed file system

Next, create a WDL, and its JSON input like so:

_Sample WDL_ 

```wdl
task sparkjob_with_yarn_cluster {
        File input_jar
        String input_1
        String output_base
        String entry_point
        Int cores
        String memory

        command {
                ${input_jar} ${input_1} ${output_base}
        }

        runtime {
                appMainClass: "${entry_point}"
                executorMemory: "${memory}"
                executorCores: "${cores}"
        }
        output {
                File out = "${output_base}"
          }
	}
```

and its accompanying JSON input as:

```json
{
	"sparkWithYarnCluster.sparkjob_with_yarn_cluster.memory": "4G",
	"sparkWithYarnCluster.sparkjob_with_yarn_cluster.output_base":"/mnt/lustre/hadoop/home/yarn_cluster_output",
	"sparkWithYarnCluster.sparkjob_with_yarn_cluster.entry_point": "com.org.spark.poc.nfs.SparkVowelLine",
	"sparkWithYarnCluster.sparkjob_with_yarn_cluster.cores": "12",
	"sparkWithYarnCluster.sparkjob_with_yarn_cluster.input_1": "/mnt/lustre/hadoop/home/inputfiles/sample.txt",
	"sparkWithYarnCluster.sparkjob_with_yarn_cluster.input_jar": "/mnt/lustre/hadoop/home/inputjars/spark_hdfs.jar"
}
```