####

# Installing the Cromwell To Use Local Scratch Device
#### These instructions are a community contribution

### In the process of being updated as of 2024-02

In this section we will install the Cromwell Workflow Management system and configure it, so it can use the local scratch device on the compute nodes.
(these installations are done in a ```centos 8``` enviornment)

### 1. In order to install Cromwell, the `sbt` build tool is required, so we will install this now.

##### a. Add the sbt online repository

```hocon
dnf config-manager --add-repo https://www.scala-sbt.org/sbt-rpm.repo
```

##### b. Install SBT build tool

```hocon
dnf -y install sbt
```

### 2. Now we can download Cromwell from git repository

##### a. Create a directory into which we will install Cromwell

```hocon
mkdir cromwell
```

##### b. Now, we set the owner of that directory to the cromwell user and group and change the permissions in order to allow access to all users.

```hocon
chmod 775 -R cromwell
chown -R cromwell:cromwell cromwell
```

##### c. Login to the cromwell user account.

```hocon
su - cromwell
```

##### d. We will now use git to clone the Cromwell git repository

```hocon
cd cromwell
git clone https://github.com/broadinstitute/cromwell.git
```

##### e. This guide was tested and validated with `version 52` of cromwell, so we will checkout that version.

```hocon
cd cromwell
git checkout 52
```

### 3. We will now configure Cromwell to use the local NVMe disk as scratch space

##### a. Open the `backend/src/main/scala/cromwell/backend/RuntimeEnvironment.scala` file for editing

##### b. Comment out line 3, so it reads:

```hocon
//import java.util.UUID
```

##### c. Remove lines 23 through 27, i.e. the following lines:

```hocon
val tempPath: String = {
val uuid = UUID.randomUUID().toString
val hash = uuid.substring(0, uuid.indexOf('-'))
callRoot.resolve(s"tmp.$hash").pathAsString
}
```

##### d. Add the following text in line 23 now:

```hocon
def tempPath: String = "/genomics_local"
```

##### e. Save and exit the file

##### f. Open the file `backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala` for editing

##### g. Go to line number 380. You will see the following content:

```hocon
|export _JAVA_OPTIONS=-Djava.io.tmpdir="$$tmpDir"
|export TMPDIR="$$tmpDir"
```


##### h. Replace those two line with the following text:

```hocon
|mkdir -p $$tmpDir/tmp.$$$$
|export _JAVA_OPTIONS=-Djava.io.tmpdir="$$tmpDir/tmp.$$$$"
|export TMPDIR="$$tmpDir/tmp.$$$$
```


##### i. Save and exit the file.You can also execute the following commands:

```hocon
sed -i "s/^\(import\ java.util.UUID\)/\/\/\1/" \
backend/src/main/scala/cromwell/backend/RuntimeEnvironment.scala
sed -i '23,27d' \
backend/src/main/scala/cromwell/backend/RuntimeEnvironment.scala
sed -i '23i \ \ \ \ \ \ def tempPath: String = \"/genomics_local\"' \
backend/src/main/scala/cromwell/backend/RuntimeEnvironment.scala
sed 's/\(\s*|export _JAVA_OPTIONS.*\)\"/\ \ \ \ \ \ \ \ |mkdir -p \$\$tmpDir\/tmp\.\$\$\$\$\n\1\"/' \
backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala
sed 's/\(\s*|export _JAVA_OPTIONS.*tmpDir\)\"/\1\/tmp\.\$\$\$\$\"/' \
backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala
sed 's/\(\s*|export TMPDIR=.*tmpDir\)\"/\1\/tmp\.\$\$\$\$\"/' \
backend/src/main/scala/cromwell/backend/standard/StandardAsyncExecutionActor.scala
```

### 4. Now we can build the patched Cromwell

```hocon
sbt clean
sbt assembly
```

### 5. When the build was successful, we can move the new jar file into the cromwell directory

```hocon
cp server/target/scala-2.13/cromwell-52-*-SNAP.jar \
cromwell/cromwell-52-fix.jar
```

## Configure the Execution Environment for Cromwell

In this section we will configure the Cromwell execution environment.We will use the default configuration file as a starting point and change it to reflect our needs.
*• Configure MariaDB as the database Cromwell should use
*• Add the `SLURM backends`, so Cromwell can schedule job using `SLURM`

### 1. First, we download the default Cromwell `reference.conf` configuration file.

```hocon
wget https://raw.githubusercontent.com/broadinstitute/cromwell/52_hotfix/core/\
src/main/resources/reference.conf cromwell/
```

### 2. We will now add `SLURM` as the backend Cromwell should use

##### a. Open the `cromwell/reference.conf` file for editing

##### b. Goto line 479 which should read:

```hocon
default = "Local"
```

##### c. Change the `"Local"` in that line to `"SLURM"`

##### d. Remove the following 5 lines, i.e. lines number 480 through 484

```hocon
providers {
Local {
actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
config {
include required(classpath("reference_local_provider_config.inc.conf"))
```

##### e. Now add the following text after line 479, i.e. after the line reading `default ="SLURM"`. Ensure that the lines that show line-breaks in this document are, in fact, single lines in `reference.conf`

```hocon
providers {
SLURM {
#Modifying temp directory to write to local disk
temporary-directory = "$(/genomics_local/)"
actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
config {
root = "cromwell-slurm-exec"
runtime-attributes = """
Int runtime_minutes = 600
Int cpu = 2
Int memory_mb = 1024
String queue = "all"
String? docker
"""
submit = """
sbatch -J ${job_name} -D ${cwd} -o ${out} -e ${err} -t ${runtime_minutes} -p ${queue} ${"-c " +
cpu} --mem ${memory_mb} --wrap "/bin/bash ${script}"
"""
kill = "scancel ${job_id}"
check-alive = "squeue -j ${job_id}"
job-id-regex = "Submitted batch job (\\d+).*"
}
}
SLURM-BWA {
temporary-directory = "$(/genomics_local/)"
actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
config {
root = "cromwell-slurm-exec"
runtime-attributes = """
Int runtime_minutes = 600
Int cpu = 2
Int memory_mb = 1024
String queue = "bwa"
String? docker
"""
submit = """
```
