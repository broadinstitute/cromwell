#!/bin/bash
#
# test_upgrade.sh
#
# What this script does:
# - Starts Cromwell using the "previous" cromwell JAR/configuration
# - Sends 3 jobs to Cromwell (2 batched, one single)
# - Waits 6 minutes
# - Sends the same 3 jobs again, in the same way
# - Shutdown Cromwell with jobs still running
# - Restarts Cromwell using the "new" cromwell JAR/configuration
# - Sends the same three jobs again
# - Waits for everything to complete
#
# What it doesn't do (yet... but maybe for the C27 release!):
# - Guarantee that any jobs are in any specific states at the point that it shuts down.
# - Check the status of jobs after completing.
# - Check that operations IDs weren't duplicated after restarting.
# - DRY it up with functions replacing the outrageous copy/pasting.
#

PREVIOUS_CROMWELL_JAR=cromwell-25-31ae549-SNAP.jar
PREVIOUS_CROMWELL_CONF=jes.conf

NEW_CROMWELL_JAR=cromwell-26-88630db-SNAP.jar
NEW_CROMWELL_CONF=jes.conf

OUT_DIR=out

PREVIOUS_CROMWELL_LOG="$OUT_DIR/previous_version.log"
NEW_CROMWELL_LOG="$OUT_DIR/new_version.log"

START_DT=$(date "+%Y-%m-%dT%H:%M:%S.000-04:00")
echo "Starting $0 at: $START_DT"
# Now we've printed it, make it URL-safe
START_DT="${START_DT//:/%3A}"

pgrep -q cromwell
ALREADY_RUNNING=$?
if [ $ALREADY_RUNNING -ne 1 ];
then
	echo "Oops! Cromwell (🐖 ) is already running!"
	pgrep cromwell
	exit 1
fi

rm -r $OUT_DIR
mkdir $OUT_DIR

echo -n "Starting Cromwell (🐖 )..."
java -Dconfig.file="$PREVIOUS_CROMWELL_CONF" -jar "$PREVIOUS_CROMWELL_JAR" server &> "$PREVIOUS_CROMWELL_LOG" &
CROMWELL_PID=$!
echo "started (PID: $CROMWELL_PID)."

echo -n "Waiting for (🐖 ) API..."
READY=-1
while [ $READY -ne 0 ];
do
	sleep 1
	curl -X GET --header "Accept: application/json" "http://localhost:8000/api/workflows/v1/backends" &>/dev/null
	READY=$?
done
echo "ready."

for i in 1 2
do
	echo -n "submitting job ${i}a and ${i}b..."
	RESULT_FILE="$OUT_DIR/submitResult_${i}_ab.json"
	curl -X POST --header "Accept: application/json" "http://localhost:8000/api/workflows/v1/batch" \
	  -F workflowSource=@scatter_files.wdl \
	  -F workflowInputs=@scatter_files_input_part1_ab.json \
	  -F workflowInputs_2=@scatter_files_input_part2.json \
	  -F workflowInputs_3=@scatter_files_input_part3.json \
	  -F workflowInputs_4=@scatter_files_input_part4.json \
	  -F workflowInputs_5=@scatter_files_input_part5.json \
	  -F workflowInputs_6=@scatter_files_input_part6.json \
	  -F workflowOptions=@defaultDocker.json \
	  -F customLabels=@custom_labels.json &> $"RESULT_FILE"
    echo "done (Response in: $RESULT_FILE)."
	echo -n "submitting job ${i}c..."
	RESULT_FILE="$OUT_DIR/submitResult_${i}_c.json"
	curl -X POST --header "Accept: application/json" "http://localhost:8000/api/workflows/v1/batch" \
	  -F workflowSource=@scatter_files.wdl \
	  -F workflowInputs=@scatter_files_input_part1_c.json \
	  -F workflowInputs_2=@scatter_files_input_part2.json \
	  -F workflowInputs_3=@scatter_files_input_part3.json \
	  -F workflowInputs_4=@scatter_files_input_part4.json \
	  -F workflowInputs_5=@scatter_files_input_part5.json \
	  -F workflowInputs_6=@scatter_files_input_part6.json \
	  -F workflowOptions=@defaultDocker.json \
	  -F customLabels=@custom_labels.json &> $"RESULT_FILE"
    echo "done (Response in: $RESULT_FILE)."
	[ "$i" -eq 1 ] && sleep 360
done

# Step two: upgrade cromwell
kill $CROMWELL_PID

pgrep -q cromwell
STILL_RUNNING=$?
while [ $STILL_RUNNING -eq 0 ];
do
	echo "Waiting for Cromwell(🐖 ) to exit..."
	sleep 1
	pgrep -q cromwell
	STILL_RUNNING=$?
done

echo -n "Starting Cromwell(🐖 )... "
java -Dconfig.file=$NEW_CROMWELL_CONF -jar $NEW_CROMWELL_JAR server &> $NEW_CROMWELL_LOG &
CROMWELL_PID=$!
echo "started (PID=$CROMWELL_PID)."

echo -n "Waiting for 🐖  API..."
READY=-1
while [ "$READY" -ne "0" ];
do
	sleep 1
	curl -X GET --header "Accept: application/json" "http://localhost:8000/api/workflows/v1/backends" &>/dev/null
	READY=$?
done
echo "ready."

i=3
echo -n "submitting job ${i}a and ${i}b..."
RESULT_FILE="$OUT_DIR/submitResult_${i}_ab.json"
curl -X POST --header "Accept: application/json" "http://localhost:8000/api/workflows/v1/batch" \
  -F workflowSource=@scatter_files.wdl \
  -F workflowInputs=@scatter_files_input_part1_ab.json \
  -F workflowInputs_2=@scatter_files_input_part2.json \
  -F workflowInputs_3=@scatter_files_input_part3.json \
  -F workflowInputs_4=@scatter_files_input_part4.json \
  -F workflowInputs_5=@scatter_files_input_part5.json \
  -F workflowInputs_6=@scatter_files_input_part6.json \
  -F workflowOptions=@defaultDocker.json \
  -F customLabels=@custom_labels.json &> $"RESULT_FILE"
echo "done (Response in: $RESULT_FILE)."
echo -n "submitting job ${i}c..."
RESULT_FILE="$OUT_DIR/submitResult_${i}_c.json"
curl -X POST --header "Accept: application/json" "http://localhost:8000/api/workflows/v1/batch" \
  -F workflowSource=@scatter_files.wdl \
  -F workflowInputs=@scatter_files_input_part1_c.json \
  -F workflowInputs_2=@scatter_files_input_part2.json \
  -F workflowInputs_3=@scatter_files_input_part3.json \
  -F workflowInputs_4=@scatter_files_input_part4.json \
  -F workflowInputs_5=@scatter_files_input_part5.json \
  -F workflowInputs_6=@scatter_files_input_part6.json \
  -F workflowOptions=@defaultDocker.json \
  -F customLabels=@custom_labels.json &> $"RESULT_FILE"
echo "done (Response in: $RESULT_FILE)."

# Step 3: Wait until everything's done:
echo -n "Waiting for the run to complete at a rate of one 🐷  per minute..."
DONE=1
while [ "$DONE" -ne "0" ];
do
	CURLED=$(curl -X GET --header "Accept: application/json" "http://localhost:8000/api/workflows/v1/query?start=${START_DT}&status=Running" 2>/dev/null)
	grep -q ": \[\]" <<< "$CURLED"
	FINISHED=$?
	grep -q '"status": "fail"' <<< "$CURLED"
	ERROR=$?
	[ $ERROR -eq 0 ] && ( echo "Error: $CURLED" )
	if [ $ERROR -eq 0 ] || [ $FINISHED -eq 0 ];
	then
		DONE=0
	else
		echo -n ".🐷 ."
		sleep 60
	fi
done
echo "...done"

kill $CROMWELL_PID

# Step 4: analyse logs to make sure things worked out:
echo "Previous version's operations IDs:"
grep operations "$PREVIOUS_CROMWELL_LOG" | sed "s/.*\ - \(.*\)/\1/g"

echo
echo "New version's operations IDs:"
grep operations "$NEW_CROMWELL_LOG" | sed "s/.*\ - \(.*\)/\1/g"
# TODO...
