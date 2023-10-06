import requests
import os
import json
import random
import string
import uuid
import time

# NOTE: POETRY prefix is outlined per documentation here: https://python-poetry.org/docs/configuration#using-environment-variables
bearer_token = os.environ['POETRY_BEARER_TOKEN']
bee_name = os.environ['POETRY_BEE_NAME']
billing_project_name = os.environ['POETRY_BILLING_PROJECT_NAME']

rawls_url = f"https://rawls.{bee_name}.bee.envs-terra.bio"
leo_url = f"https://leonardo.{bee_name}.bee.envs-terra.bio"

def output_message(msg):
    current_time = time.strftime("%H:%M:%S", time.localtime())
    print(f'{current_time} - {msg}')

def handle_failed_request(response, msg, status_code=200):
    if(response.status_code != status_code):
        raise Exception(f'{response.status_code} - {msg}')

def create_workspace():
   rawls_api_call = f"{rawls_url}/api/workspaces"
   request_body= {
      "namespace": billing_project_name, # Billing project name
      "name": f"api-workspace-{''.join(random.choices(string.ascii_lowercase, k=5))}", # workspace name
      "attributes": {}}
   
   create_workspace_response = requests.post(url=rawls_api_call, 
                                             json=request_body, 
                                             headers={"Authorization": f"Bearer {bearer_token}"}
   ).json()

   create_workspace_data = json.loads(json.dumps(create_workspace_response))
   workspace_id = create_workspace_data['workspaceId']

   output_message(f"Enabling CBAS for workspace {workspace_id}")
   activate_cbas_request = f"{leo_url}/api/apps/v2/{workspace_id}/terra-app-{str(uuid.uuid4())}"
   cbas_request_body = {
      "appType": "CROMWELL"
   } 
        
   response = requests.post(url=activate_cbas_request, json=cbas_request_body, 
                            headers={"Authorization": f"Bearer {bearer_token}"})
   # will return 202 or error
   handle_failed_request(response, "Error activating CBAS", 202)
   return create_workspace_data

# GET WDS OR CROMWELL ENDPOINT URL FROM LEO
def get_app_url(workspaceId, app):
    uri = f"{leo_url}/api/apps/v2/{workspaceId}?includeDeleted=false"

    headers = {"Authorization": bearer_token,
               "accept": "application/json"}

    response = requests.get(uri, headers=headers)
    status_code = response.status_code

    if status_code != 200:
        return response.text
    output_message("Successfully retrieved details.")
    response = response.json()

    app_url = ""
    app_type = "CROMWELL" if app != 'wds' else app.upper()
    output_message(f"App type: {app_type}")
    for entries in response: 
        if entries['appType'] == app_type and entries['proxyUrls'][app] is not None:
            print(entries['status'])
            if(entries['status'] == "PROVISIONING"):
                output_message(f"{app} is still provisioning")
                break
            output_message(f"App status: {entries['status']}")
            app_url = entries['proxyUrls'][app]
            break 

    if app_url is None: 
        output_message(f"{app} is missing in current workspace")
    else:
        output_message(f"{app} url: {app_url}")

    return app_url

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
    success_statuses = ['Succeeded']
    throw_exception_statuses = ['Aborted', 'Failed'] # Are there other statuses that should throw an exception?
    
    current_running_workflow_count = 0
    while workflow_ids:
        if max_retries == 0:
            raise Exception(f"Workflow(s) did not finish running within retry window ({max_retries} retries)")
        
        workflow_id = workflow_ids.pop()
        workflow_metadata = get_workflow_information(app_url, workflow_id)
        workflow_status = workflow_metadata['status']

        if(workflow_status in throw_exception_statuses):
            raise Exception(f"Exception raised: Workflow {workflow_id} reporting {workflow_status} status")
        if workflow_status in success_statuses:
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

    output_message("Workflow(s) submission and completion successful")

def deleteApps(workspace_id):
    delete_url = f"{leo_url}/api/apps/v2/{workspace_id}/deleteAll"
    response = requests.delete(url=delete_url, 
                                headers={"Authorization": f"Bearer {bearer_token}"})
    handle_failed_request(response, f"Error deleting apps for workspace {workspace_id}", 202)
    output_message(f"Apps successfully deleted for workspace {workspace_id}")

def deleteWorkspace(workspace_namespace, workspace_name):
    if workspace_namespace and workspace_name:
        rawls_api_call = f"{rawls_url}/api/workspaces/{workspace_namespace}/{workspace_name}"
        response = requests.delete(url=rawls_api_call, 
                                headers={"Authorization": f"Bearer {bearer_token}"}
        )
        handle_failed_request(response, f"Error deleting workspace {workspace_namespace} - {workspace_name}", 204)
        output_message(f"Workspace {workspace_name} successfully deleted")

def start():
    workspace_namespace = ""
    workspace_name = ""
    # Giving workflow 3 minutes to complete
    sleep_timer = 60 * 3
    try:
        created_workspace = create_workspace()
        workspace_id = created_workspace['workspaceId']
        workspace_namespace = created_workspace['namespace']
        workspace_name = created_workspace['name']
        time.sleep(60 * 20) # Added an sleep here to give the workspace time to provision
        app_url = get_app_url(workspace_id, 'cromwell')

        # This chunk of code only executes one workflow
        # Would like to modify this down the road to execute and store references for multiple workflows
        workflow_response = submit_workflow_to_cromwell(app_url, "Run Workflow Test")
        #Will need to update this when swapping out hello wdl with fetch_sra_to_bam (20 min?)
        output_message(f'Executing sleep for {sleep_timer} seconds to allow workflow(s) to complete')
        time.sleep(sleep_timer)

        # This chunk of code supports checking one or more workflows
        # Probably won't require too much modification if we want to run additional submission tests
        workflow_ids = [workflow_response['id']]
        get_completed_workflow(app_url, workflow_ids)
    except Exception as e:
        raise e
    finally:
        deleteApps(workspace_id)
        time.sleep(60 * 3) # Not sure if this is necessary
        deleteWorkspace(workspace_namespace, workspace_name)
        time.sleep(60 * 3) # Not sure if this is necessary
        output_message("Workspace cleanup complete")