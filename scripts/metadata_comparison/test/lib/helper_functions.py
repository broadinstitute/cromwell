#!/usr/bin/env python3

def read_resource(filename):
    path = f'test/resources/{filename}'
    with open(path, 'r') as file:
        data = file.read()
    return data