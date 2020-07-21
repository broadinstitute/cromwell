#!/usr/bin/env python3
#
# Initializer for various PAPI clients which only creates clients for specific APIs when we need them
#

from googleapiclient.discovery import build as google_client_build, Resource
from google.auth.credentials import Credentials
from metadata_comparison.lib.operation_ids import operation_id_to_api_version
from typing import Mapping, Any
import logging

logger = logging.getLogger('metadata_comparison.lib.papi.PapiClients')


class PapiClients:
    clients = {}

    def __init__(self, credentials: Credentials) -> None:
        self.credentials = credentials

    def __get_client(self, api_version: str) -> Resource:
        """Gets the relevant client for accessing a PAPI API, or makes a new instance if necessary"""
        if api_version not in self.clients:
            self.clients[api_version] = self.__make_client(api_version)
        return self.clients[api_version]

    def __make_client(self, api_version: str) -> Resource:
        """Makes a new client for accessing a specified PAPI API"""
        if api_version in ['v1alpha2', 'v2alpha1']:
            return google_client_build('genomics', api_version, credentials=self.credentials)
        elif api_version == 'v2beta':
            return google_client_build('lifesciences', api_version, credentials=self.credentials)
        else:
            raise Exception(f'Unsupported client api_version: "{api_version}"')

    @staticmethod
    def __read_papi_v1_operation_metadata(operation_id: str, genomics_v1_client: Resource) -> Mapping[str, Any]:
        """Reads the operations metadata for a pipelines API v1 job ID. Returns a python dict"""
        logger.info(f'Reading PAPI v1 operation metadata for {operation_id}...')
        result = genomics_v1_client.operations().get(name=operation_id).execute()
        return result

    @staticmethod
    def __read_papi_v2alpha1_operation_metadata(operation_id: str, genomics_v2alpha1_client: Resource) -> Mapping[str, Any]:
        """Reads the operations metadata for a pipelines API v2alpha1 job ID. Returns a python dict"""
        logger.info(f'Reading PAPI v2alpha1 operation metadata for {operation_id}...')
        result = genomics_v2alpha1_client.projects().operations().get(name=operation_id).execute()
        return result

    @staticmethod
    def __read_papi_v2beta_operation_metadata(operation_id: str, genomics_v2beta_client: Resource) -> Mapping[str, Any]:
        """Reads the operations metadata for a pipelines API v2beta job ID. Returns a python dict"""
        logger.info(f'Reading PAPI v2beta operation metadata for {operation_id}...')
        result = genomics_v2beta_client.projects().locations().operations().get(name=operation_id).execute()
        return result

    def request_operation_metadata(self, operation_id: str) -> Mapping[str, Any]:
        """
        Reads the operations metadata for any supported pipelines API version.
        Returns a python dict
        """
        api_version = operation_id_to_api_version(operation_id)
        client = self.__get_client(api_version)
        if api_version == 'v1alpha2':
            return self.__read_papi_v1_operation_metadata(operation_id, client)
        elif api_version == 'v2alpha1':
            return self.__read_papi_v2alpha1_operation_metadata(operation_id, client)
        elif api_version == 'v2beta':
            return self.__read_papi_v2beta_operation_metadata(operation_id, client)
        else:
            raise Exception(f'Unsupported client api_version in request_operation_metadata: "{api_version}"')
