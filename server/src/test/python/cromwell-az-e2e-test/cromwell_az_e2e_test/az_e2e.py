import requests
import os
import json
import random
import string
import uuid
import time

# NOTE: defining env variables for appUrl for now, will remove when provisioning ticket is in progress
app_url = os.environ['POETRY_APP_URL']
bearer_token = os.environ['POETRY_BEARER_TOKEN']

# Commenting out GHS specific variables for now
# bee_name = os.environ['BEE_NAME']
# billing_project_name = os.environ['BILLING_PROJECT_NAME']
# number_of_workspaces = 1
# wds_upload=False
# cbas_submit_workflow=False
# number_of_workflows_to_kick_off = 1

# rawls_url = f"https://rawls.{bee_name}.bee.envs-terra.bio"
# leo_url = f"https://leonardo.{bee_name}.bee.envs-terra.bio"

def output_message(msg):
    current_time = time.strftime("%H:%M:%S", time.localtime())
    print(f'{current_time} - {msg}')

def handle_failed_request(response, msg, status_code=200):
    if(response.status_code != status_code):
        raise Exception(f'{response.status_code} - {msg}')

# Commenting out workspace/provisioning helpers for now
# def create_workspace():
#    rawls_api_call = f"{rawls_url}/api/workspaces"
#    request_body= {
#       "namespace": billing_project_name, # Billing project name
#       "name": f"api-workspace-{''.join(random.choices(string.ascii_lowercase, k=5))}", # workspace name
#       "attributes": {}}
   
#    create_workspace_response = requests.post(url=rawls_api_call, 
#                                              json=request_body, 
#                                              headers={"Authorization": f"Bearer {bearer_token}"}
#    ).json()

#    create_workspace_data = json.loads(json.dumps(create_workspace_response))
#    workspaceId = create_workspace_data['workspaceId']

#    print(f"Enabling CBAS for workspace {workspaceId}")
#    activate_cbas_request = f"{leo_url}/api/apps/v2/{workspaceId}/terra-app-{str(uuid.uuid4())}"
#    cbas_request_body = {
#       "appType": "CROMWELL"
#    } 
        
#    response = requests.post(url=activate_cbas_request, json=cbas_request_body, 
#                             headers={"Authorization": f"Bearer {bearer_token}"})
#    # will return 202 or error
#    handle_failed_request(response, "Error activating CBAS", 202)
   
#    print(response)
#    return workspaceId

# # GET WDS OR CROMWELL ENDPOINT URL FROM LEO
# def get_app_url(workspaceId, app):
#     """"Get url for wds/cbas."""
#     uri = f"{leo_url}/api/apps/v2/{workspaceId}?includeDeleted=false"

#     headers = {"Authorization": bearer_token,
#                "accept": "application/json"}

#     response = requests.get(uri, headers=headers)
#     status_code = response.status_code

#     if status_code != 200:
#         return response.text
#     print("Successfully retrieved details.")
#     response = response.json()

#     app_url = ""
#     app_type = "CROMWELL" if app != 'wds' else app.upper()
#     print(f"App type: {app_type}")
#     for entries in response: 
#         if entries['appType'] == app_type and entries['proxyUrls'][app] is not None:
#             print(entries['status'])
#             if(entries['status'] == "PROVISIONING"):
#                 print(f"{app} is still provisioning")
#                 break
#             print(f"App status: {entries['status']}")
#             app_url = entries['proxyUrls'][app]
#             break 

#     if app_url is None: 
#         print(f"{app} is missing in current workspace")
#     else:
#         print(f"{app} url: {app_url}")

#     return app_url

def submit_workflow_to_cromwell(app_url, workflow_test_name):
    absolute_file_path = os.path.dirname(__file__)
    workflow_source_path = os.path.join(absolute_file_path, 'workflow_files/hello.wdl')
    workflow_inputs_path = os.path.join(absolute_file_path, 'workflow_files/hello.inputs')
    workflow_endpoint = f'{app_url}/api/workflows/v1'
    headers = {"Authorization": f'Bearer {bearer_token}',
              "accept": "application/json",
    }
    files = {
        'workflowSource': ('hello.wdl', open(workflow_source_path, 'rb'), 'application/octet-stream'),
        'workflowInputs': ('hello.inputs', open(workflow_inputs_path, 'rb'), 'application/octet-stream'),
        'workflowType': 'WDL',
        'workflowTypeVersion': '1.0'
    }
    response = requests.post(workflow_endpoint, headers=headers, files=files)
    handle_failed_request(response, f"Error submitting workflow to Cromwell for {workflow_test_name}", 201)
    output_message(response.json())
    return response.json()

def get_workflow_information(app_url, workflow_id):
    workflow_endpoint = f'{app_url}/api/workflows/v1/{workflow_id}/metadata'
    headers = {"Authorization": f'Bearer {bearer_token}',
              "accept": "application/json"}
    output_message(f"Fetching workflow metadata for {workflow_id}")
    response = requests.get(workflow_endpoint, headers=headers)
    handle_failed_request(response, f"Error fetching workflow metadata for {workflow_id}")
    output_message(response.json())
    return response.json()

def get_completed_workflow(app_url, workflow_ids, max_retries=4):
    target_statuses = ['Succeeded', 'Failed']
    current_running_workflow_count = 0
    while workflow_ids:
        if max_retries == 0:
            raise Exception(f"Workflow(s) did not finish running within retry window ({max_retries} retries)")
        workflow_id = workflow_ids.pop()
        workflow_metadata = get_workflow_information(app_url, workflow_id)
        if workflow_metadata['status'] in target_statuses:
            output_message(f"{workflow_id} finished running. Status: {workflow_metadata['status']}")
        else:
            workflow_ids.append(workflow_id)
            current_running_workflow_count += 1
        if current_running_workflow_count == len(workflow_ids):
            if current_running_workflow_count == 0:
                output_message("Workflow(s) finished running")
            else:
                # Reset current count to 0 for next retry
                # Decrement max_retries by 1
                # Wait 2 minutes before checking workflow statuses again (adjust as needed)
                output_message(f"These workflows have yet to return a completed status: [{', '.join(workflow_ids)}]")
                max_retries -= 1
                current_running_workflow_count = 0
                time.sleep(60 * 2)

def start():
    sleep_timer = 60 * 3

    # Commenting out the workspace/provisioning steps for now
    # workspace_id = create_workspace()
    # time.sleep(60 * 20) # Added an sleep here to give the workspace time to provision
    # app_url = get_app_url(workspace_id, 'cromwell')

    # This chunk of code only executes one workflow
    # Would like to modify this down the road to execute and store references for multiple workflows
    workflow_response = submit_workflow_to_cromwell(app_url, "Run Workflow Test")
    #Giving workflow 10 minutes to complete
    #Will need to update this when swapping out hello wdl with fetch_sra_to_bam (20 min?)
    output_message(f'Executing sleep for {sleep_timer} seconds to allow workflow(s) to complete')
    time.sleep(sleep_timer)

    # This chunk of code supports checking one or more workflows
    # Probably won't require too much modification if we want to run additional submission tests
    workflow_ids = [workflow_response['id']]
    get_completed_workflow(app_url, workflow_ids)
    output_message("Workflow(s) submission and completion successful")