#!/usr/bin/env bash

set -e

removeCromwellJar() {
    gsutil rm "${JAR_GCS_PATH}" || true
}

printTravisHeartbeat() {
    # Sleep one minute between printouts, but don't zombie for more than two hours
    for ((i=0; i < 120; i++)); do
        sleep 60
        printf "…"
    done &
    TRAVIS_HEARTBEAT_PID=$!
}

killTravisHeartbeat() {
    if [ -n "${TRAVIS_HEARTBEAT_PID+set}" ]; then
        kill ${TRAVIS_HEARTBEAT_PID} || true
    fi
}

exitScript() {
    killTravisHeartbeat
    removeCromwellJar
}

trap exitScript EXIT
printTravisHeartbeat

set -x

# Unpack our credentials and such
openssl aes-256-cbc -K "$encrypted_5b9e82629fa8_key" -iv "$encrypted_5b9e82629fa8_iv" -in src/bin/travis/resources/jesConf.tar.enc -out jesConf.tar -d
tar xvf jesConf.tar

# Do a bunch of crap to enable gsutil. It's possible this is overkill but it doesn't take long anyways
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 1397BC53640DB551
CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)"
echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | sudo tee /etc/apt/sources.list.d/google-cloud-sdk.list
curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo apt-key add -
sudo apt-get install google-cloud-sdk
export PYTHONPATH="/usr/lib/python2.7/site-packages:/usr/local/lib/python2.7/site-packages:/usr/local/lib/python2.7/dist-packages"
export CONFIGURE_OPTS="--enable-unicode=ucs4"
pyenv install 2.7.10
pyenv global 2.7.10
sudo pip install --upgrade pip
sudo pip install pyopenssl ndg-httpsclient pyasn1 --upgrade
export CLOUDSDK_PYTHON_SITEPACKAGES=1
gcloud auth activate-service-account --key-file=broad-dsde-cromwell-dev-d9c443bb4a94.json "$CROMWELL_SERVICE_ACCOUNT"

echo "RUNNING TRAVIS CENTAUR"
sbt assembly

# Update the inputs file with stuff specific to this run
sed -i "s/CENTAUR_BRANCH/${CENTAUR_BRANCH}/g" src/bin/travis/resources/centaur.inputs
CROMWELL_JAR=cromwell_${TRAVIS_BUILD_ID}.jar
sed -i "s/CROMWELL_JAR/${CROMWELL_JAR}/g" src/bin/travis/resources/centaur.inputs

# Upload the built Cromwell jar to GCS so we can use it in our centaur test. Set an exit trap to clean it up on failure
JAR_GCS_PATH=gs://cloud-cromwell-dev/travis-centaur/${CROMWELL_JAR}
gsutil cp target/scala-2.11/cromwell-*.jar "${JAR_GCS_PATH}"

java -Dconfig.file=./jes.conf -jar target/scala-2.11/cromwell-*.jar run src/bin/travis/resources/centaur.wdl src/bin/travis/resources/centaur.inputs | tee log.txt
EXIT_CODE="${PIPESTATUS[0]}"

# The perl code below is to remove our lovely color highlighting
export WORKFLOW_ID=`grep "SingleWorkflowRunnerActor: Workflow submitted " log.txt | perl -pe 's/\e\[?.*?[\@-~]//g' | cut -f7 -d" "`

# Grab the Centaur log from GCS and cat it so we see it in the main travis log. 
export CENTAUR_LOG_PATH="gs://cloud-cromwell-dev/cromwell_execution/travis/centaur_workflow/${WORKFLOW_ID}/call-centaur/cromwell_root/logs/centaur.log"
gsutil cp ${CENTAUR_LOG_PATH} centaur.log
cat centaur.log
echo "More logs for this run are available at https://console.cloud.google.com/storage/browser/cloud-cromwell-dev/cromwell_execution/travis/centaur_workflow/${WORKFLOW_ID}/call-centaur/"

exit "${EXIT_CODE}"
