#!/usr/bin/env python3
#
# Lazy initializer for various PAPI clients so that we only create them when we need them
#

from googleapiclient.discovery import build as google_client_build
from metadata_comparison.lib.operation_ids import *
import logging

logger = logging.getLogger('metadata_comparison.lib.papi.PapiClients')

class PapiClients:
    clients = {}

    def __init__(self, credentials):
        self.credentials = credentials


    def get_client(self, api_version):
        if api_version not in self.clients:
            self.clients[api_version] = self.make_client(api_version)
        return self.clients[api_version]


    def make_client(self, api_version):
        if api_version in [ 'v1alpha2', 'v2alpha1']:
            return google_client_build('genomics', api_version, credentials = self.credentials)
        elif api_version == 'v2beta':
            return google_client_build('lifesciences', api_version, credentials = self.credentials)
        else:
            raise Exception(f'Unsupported client api_version: "{api_version}"')


    def read_papi_v1_operation_metadata(self, operation_id, genomics_v1_client):
        """Reads the operations metadata for a pipelines API v1 job ID. Returns a python dict"""
        logger.info(f'Reading PAPI v1 operation metadata for {operation_id}...')
        result = genomics_v1_client.operations().get(name=operation_id).execute()
        return result


    def read_papi_v2alpha1_operation_metadata(self, operation_id, genomics_v2alpha1_client):
        """Reads the operations metadata for a pipelines API v2alpha1 job ID. Returns a python dict"""
        logger.info(f'Reading PAPI v2alpha1 operation metadata for {operation_id}...')
        result = genomics_v2alpha1_client.projects().operations().get(name=operation_id).execute()
        return result


    def read_papi_v2beta_operation_metadata(self, operation_id, genomics_v2beta_client):
        """Reads the operations metadata for a pipelines API v2beta job ID. Returns a python dict"""
        logger.info(f'Reading PAPI v2beta operation metadata for {operation_id}...')
        result = genomics_v2beta_client.projects().locations().operations().get(name=operation_id).execute()
        return result


    def request_operation_metadata(self, operation_id):
        """
        Reads the operations metadata for any supported pipelines API version.
        Returns a python dict
        """
        api_version = operation_id_to_api_version(operation_id)
        client = self.get_client(api_version)
        if api_version == 'v1alpha2':
            return self.read_papi_v1_operation_metadata(operation_id, client)
        elif api_version == 'v2alpha1':
            return self.read_papi_v2alpha1_operation_metadata(operation_id, client)
        elif api_version == 'v2beta':
            return self.read_papi_v2beta_operation_metadata(operation_id, client)
        else:
            raise Exception(f'Unsupported client api_version in request_operation_metadata: "{api_version}"')
