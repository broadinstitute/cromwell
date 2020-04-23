#!/usr/bin/env python3
#
# Lazy initializer for various PAPI clients so that we only create them when we need them
#

from googleapiclient.discovery import build as google_client_build

class PapiClients:
    def __init__(self, credentials):
        self.credentials = credentials

    clients = {}

    def get_client(self, api_version):
        if api_version not in self.clients:
            self.clients[api_version] = self.make_client(api_version)
        return self.clients[api_version]


    def make_client(self, api_version):
        if 'api_version' in [ 'v1alpha2', 'v2alpha1']:
            return google_client_build('genomics', 'v2alpha1', credentials = credentials)
        else raise Exception(f'Unsupported client api_version: "{api_version}"')
