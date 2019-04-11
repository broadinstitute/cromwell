{
    "cwlVersion": "v1.0", 
    "$graph": [
        {
            "inputs": [
                {
                    "type": "string", 
                    "id": "#main/bam"
                }, 
                {
                    "type": "#capture_kit.yml/capture_kit", 
                    "id": "#main/capture_kit"
                }
            ], 
            "requirements": [
                {
                    "class": "SchemaDefRequirement", 
                    "types": [
                        {
                            "fields": [
                                {
                                    "type": "string", 
                                    "name": "#capture_kit.yml/capture_kit/bait"
                                }
                            ], 
                            "type": "record", 
                            "name": "#capture_kit.yml/capture_kit"
                        }
                    ]
                }
            ], 
            "outputs": [
                {
                    "outputSource": "#main/touch_bam/output", 
                    "type": "File", 
                    "id": "#main/output_bam"
                }
            ], 
            "class": "Workflow", 
            "steps": [
                {
                    "out": [
                        {
                            "id": "#main/touch_bam/output"
                        }
                    ], 
                    "run": "#touch.cwl", 
                    "id": "#main/touch_bam", 
                    "in": [
                        {
                            "source": "#main/bam", 
                            "id": "#main/touch_bam/input"
                        }
                    ]
                }
            ], 
            "id": "#main"
        }, 
        {
            "inputs": [
                {
                    "inputBinding": {
                        "position": 0
                    }, 
                    "type": "string", 
                    "id": "#touch.cwl/input"
                }
            ], 
            "requirements": [
                {
                    "dockerPull": "ubuntu:bionic-20180426", 
                    "class": "DockerRequirement"
                }
            ], 
            "outputs": [
                {
                    "outputBinding": {
                        "glob": "$(inputs.input)"
                    }, 
                    "type": "File", 
                    "id": "#touch.cwl/output"
                }
            ], 
            "baseCommand": [
                "touch"
            ], 
            "class": "CommandLineTool", 
            "id": "#touch.cwl"
        }
    ]
}
