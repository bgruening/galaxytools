{
    "$graph": [
        {
            "class": "CommandLineTool",
            "inputs": [
                {
                    "type": "string",
                    "inputBinding": {
                        "position": 1,
                        "valueFrom": "A_$(inputs.in1)_B_$(inputs.in1)_C_$(inputs.in1)"
                    },
                    "id": "#tool1.cwl/in1"
                }
            ],
            "baseCommand": "echo",
            "arguments": [
                {
                    "valueFrom": "$(runtime)"
                }
            ],
            "stdout": "out.txt",
            "requirements": [
                {
                    "expressionLib": [
                        "var foo = function(x) {\n    return 2 * x\n}\n\nvar bar = function(n, x) {\n    return `{n} engineers walk into a {x}`\n}"
                    ],
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "id": "#tool1.cwl",
            "outputs": [
                {
                    "type": "string",
                    "outputBinding": {
                        "glob": "out.txt",
                        "loadContents": true,
                        "outputEval": "$(self)_D_$(runtime)"
                    },
                    "id": "#tool1.cwl/out1"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "inputs": [
                {
                    "type": "#testtypes.yml/my_boolean_array",
                    "inputBinding": {
                        "position": 1,
                        "valueFrom": "A_$(inputs.in1)_B_$(inputs.in1)_C_$(inputs.in1)"
                    },
                    "id": "#tool2.cwl/in1"
                }
            ],
            "baseCommand": "echo",
            "arguments": [
                {
                    "valueFrom": "$(runtime)"
                }
            ],
            "outputs": [
                {
                    "type": "string",
                    "outputBinding": {
                        "glob": "out.txt",
                        "loadContents": true,
                        "outputEval": "$(self)_D_$(runtime)"
                    },
                    "id": "#tool2.cwl/out1"
                }
            ],
            "stdout": "out.txt",
            "requirements": [
                {
                    "types": [
                        {
                            "name": "#testtypes.yml/my_boolean_array",
                            "type": "array",
                            "items": "boolean",
                            "label": "A boolean array"
                        },
                        {
                            "name": "#testtypes.yml/my_enum",
                            "type": "enum",
                            "symbols": [
                                "#testtypes.yml/my_enum/a",
                                "#testtypes.yml/my_enum/b",
                                "#testtypes.yml/my_enum/c"
                            ],
                            "label": "A required enum"
                        }
                    ],
                    "class": "SchemaDefRequirement"
                }
            ],
            "id": "#tool2.cwl"
        },
        {
            "class": "Workflow",
            "inputs": [
                {
                    "id": "#main/in1",
                    "type": "#testtypes.yml/my_boolean_array"
                }
            ],
            "steps": [
                {
                    "run": "#tool2.cwl",
                    "in": [
                        {
                            "source": "#main/in1",
                            "id": "#main/s1/in1"
                        }
                    ],
                    "out": [
                        "#main/s1/out1"
                    ],
                    "id": "#main/s1"
                },
                {
                    "run": "#tool1.cwl",
                    "in": [
                        {
                            "source": "#main/s1/out1",
                            "id": "#main/s2/in1"
                        }
                    ],
                    "out": [
                        "#main/s2/out1"
                    ],
                    "id": "#main/s2"
                }
            ],
            "outputs": [
                {
                    "id": "#main/out1",
                    "type": "string",
                    "outputSource": "#main/s2/out1"
                }
            ],
            "requirements": [
                {
                    "types": [
                        {
                            "$import": "#testtypes.yml/my_boolean_array"
                        },
                        {
                            "$import": "#testtypes.yml/my_enum"
                        }
                    ],
                    "class": "SchemaDefRequirement"
                }
            ],
            "id": "#main"
        }
    ],
    "cwlVersion": "v1.2"
}
