# operation_inquisitor.py
#
# Purpose: Read operations IDs from the first column of a CSV file input, and
#          output machine type and duration details for each operation (sorted
#          by project).
#
# Usage: python3 operation_inquisitor.py <<operations_list.csv>> [<<project1>> <<project 2>> ...  <<project n>>]
#
# Python Prereqs (at least, the ones which I needed to manually install... YMMV):
#   * pip3 install python-dateutil
#   * pip3 install yaml
#
#  Other Prereqs:
#   * You must have gcloud installed and configured so that it can make operations
#     describe requests for all of the operation IDs in the CSV file.

import csv
import re
import sys
import subprocess
import json
from dateutil.parser import parse as parsedate

def getOperationsMetadata(operationId):
    print('Getting operation metadata for: ' + operationId)
    completedProcess = subprocess.run(["gcloud", "alpha", "genomics", "operations", "describe", operationId, "--format", "json"],
                                      capture_output=True,
                                      encoding='ascii')
    if completedProcess.returncode == 0:
        try:
            operation_metadata = json.loads(completedProcess.stdout)
            metadata_field = operation_metadata['metadata']
            start = parsedate(metadata_field['startTime'])
            end = parsedate(metadata_field['endTime'])
            duration = (end - start).total_seconds()
            events = metadata_field['events']
            for event in events:
                if event['description'].startswith('Worker "google-pipelines-worker-'):
                    machine_type = event['details']['machineType']
                    instance = event['details']['instance']
                    zone = event['details']['zone']
            disks = metadata_field['pipeline']['resources']['virtualMachine']['disks']
            disk_details = [ disk['name'] + ' ' + str(disk['sizeGb']) + 'GB ' + disk['type'] for disk in disks ]
            return [operationId, machine_type, instance, zone, str(duration), ", ".join(disk_details)]
        except:
            print('Unable to process ' + operationId + '. Error: ' + str(sys.exc_info()[1]) + ': ' + str(sys.exc_info()[1]))
            return [operationId, "???", "???", "???", "???", "???"]
    else:
        print('Unable to process ' + operationId + '. Exit code: ' + str(completedProcess.returncode) + '. Stderr: ' + completedProcess.stderr)
        return [operationId, "???", "???", "???", "???", "???"]


# Read input CSV file name from command line:
if len(sys.argv) < 2:
    print('Usage: $ python3 ' + sys.argv[0] + ' <<operations_list.csv>> [<<project1>> <<project 2>> ...  <<project n>>]')
    sys.exit(1)
csv_in = sys.argv[1]
requested_projects = sys.argv[2:]
if len(requested_projects) == 0:
    print('Examining all rows in input CSV')
else:
    print('Filtering for the following projects in input CSV: ' + str(requested_projects))

# Read in the project IDs:
projects = {}
with open(csv_in) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        if len(row) == 0:
            break
        operationId = row[0]
        project_match = re.search('projects/([^/]*)/.*', operationId)
        if project_match is None:
            print('Skipping row ' + ', '.join(row) + ': ' + row[0] + ' does not match operation ID regex.')
        else:
            project = project_match.group(1)
            if len(requested_projects) == 0 or project in requested_projects:
                operations = projects.get(project, [])
                projects[project] = operations + [operationId]

for project in projects:
    operations = projects[project]
    print('Processing project ' + project + '. Total operations in this project: ' + str(len(operations)))
    operation_details = [ getOperationsMetadata(operation) for operation in operations ]

    output_file = 'duration_details_' + project + '.csv'
    with open(output_file, 'w+') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writerow(["Project", "OperationId", "MachineType", "Instance", "Zone", "DurationSeconds", "Disks"])
        for operation_details_entry in operation_details:
            row = [project] + operation_details_entry
            csv_writer.writerow(row)
    print('*** Duration details for project ' + project + ' written to: ' + output_file)




