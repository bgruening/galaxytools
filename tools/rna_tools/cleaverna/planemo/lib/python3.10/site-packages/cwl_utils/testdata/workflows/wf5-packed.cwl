{
    "$graph": [
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "listing": [
                        {
                            "entryname": "inputs.txt",
                            "entry": "$(inputs.in1.file.path)\n$(inputs.in1.meta.species)\n"
                        }
                    ],
                    "class": "InitialWorkDirRequirement"
                },
                {
                    "types": [
                        {
                            "name": "#recursive.yml/sample_meta",
                            "type": "record",
                            "fields": [
                                {
                                    "name": "#recursive.yml/sample_meta/sample",
                                    "type": [
                                        "null",
                                        "string"
                                    ]
                                },
                                {
                                    "name": "#recursive.yml/sample_meta/species",
                                    "type": "string"
                                }
                            ]
                        },
                        {
                            "name": "#recursive.yml/file_with_sample_meta",
                            "type": "record",
                            "fields": [
                                {
                                    "name": "#recursive.yml/file_with_sample_meta/file",
                                    "type": "File"
                                },
                                {
                                    "name": "#recursive.yml/file_with_sample_meta/meta",
                                    "type": "#recursive.yml/sample_meta"
                                }
                            ]
                        },
                        {
                            "name": "#recursive.yml/info_with_sample_meta",
                            "type": "record",
                            "fields": [
                                {
                                    "type": "string",
                                    "name": "#recursive.yml/info_with_sample_meta/comment"
                                },
                                {
                                    "type": "#recursive.yml/sample_meta",
                                    "name": "#recursive.yml/info_with_sample_meta/meta"
                                }
                            ]
                        },
                        {
                            "name": "#array.yml/sample_meta2",
                            "type": "record",
                            "fields": [
                                {
                                    "name": "#array.yml/sample_meta2/prop",
                                    "type": "string"
                                }
                            ]
                        },
                        {
                            "name": "#array.yml/study_meta",
                            "type": "array",
                            "items": "#array.yml/sample_meta2"
                        },
                        {
                            "name": "#array.yml/study_meta_too",
                            "type": "record",
                            "fields": [
                                {
                                    "type": "#array.yml/sample_meta2",
                                    "name": "#array.yml/study_meta_too/meta1"
                                },
                                {
                                    "type": "#array.yml/study_meta",
                                    "name": "#array.yml/study_meta_too/meta2"
                                }
                            ]
                        },
                        {
                            "name": "#singletype.yml/simple_record",
                            "type": "record",
                            "fields": [
                                {
                                    "type": "string",
                                    "name": "#singletype.yml/simple_record/prop"
                                }
                            ]
                        }
                    ],
                    "class": "SchemaDefRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "#recursive.yml/file_with_sample_meta",
                    "id": "#link-to-clt1.cwl/in1"
                },
                {
                    "type": "#array.yml/study_meta_too",
                    "id": "#link-to-clt1.cwl/in2"
                },
                {
                    "type": "#singletype.yml/simple_record",
                    "id": "#link-to-clt1.cwl/in3"
                },
                {
                    "type": [
                        "string",
                        "#recursive.yml/sample_meta"
                    ],
                    "id": "#link-to-clt1.cwl/in4"
                },
                {
                    "type": [
                        "null",
                        "Any"
                    ],
                    "id": "#link-to-clt1.cwl/in5"
                }
            ],
            "baseCommand": [
                "echo"
            ],
            "arguments": [
                "hello world"
            ],
            "id": "#link-to-clt1.cwl",
            "stdout": "2cd5f434d33dce1a50ec686c741fba97b41d2544",
            "hints": [
                {
                    "class": "LoadListingRequirement",
                    "loadListing": "deep_listing"
                },
                {
                    "class": "NetworkAccess",
                    "networkAccess": true
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.txt"
                    },
                    "id": "#link-to-clt1.cwl/out1"
                },
                {
                    "type": "#array.yml/study_meta_too",
                    "outputBinding": {
                        "outputEval": "$(inputs.in2)"
                    },
                    "id": "#link-to-clt1.cwl/out2"
                },
                {
                    "type": "File",
                    "id": "#link-to-clt1.cwl/out3",
                    "outputBinding": {
                        "glob": "2cd5f434d33dce1a50ec686c741fba97b41d2544"
                    }
                }
            ]
        },
        {
            "class": "Workflow",
            "inputs": [
                {
                    "type": "#recursive.yml/file_with_sample_meta",
                    "id": "#main/in1"
                },
                {
                    "type": "#array.yml/study_meta_too",
                    "id": "#main/in2"
                },
                {
                    "type": "#singletype.yml/simple_record",
                    "id": "#main/in3"
                },
                {
                    "type": [
                        "string",
                        "#recursive.yml/sample_meta"
                    ],
                    "id": "#main/in4"
                }
            ],
            "steps": [
                {
                    "run": "#link-to-clt1.cwl",
                    "in": [
                        {
                            "source": "#main/in1",
                            "id": "#main/s1/in1"
                        },
                        {
                            "source": "#main/in2",
                            "id": "#main/s1/in2"
                        },
                        {
                            "source": "#main/in3",
                            "id": "#main/s1/in3"
                        },
                        {
                            "source": "#main/in4",
                            "id": "#main/s1/in4"
                        }
                    ],
                    "out": [
                        "#main/s1/out2"
                    ],
                    "id": "#main/s1"
                }
            ],
            "outputs": [
                {
                    "id": "#main/out1",
                    "type": "#array.yml/study_meta_too",
                    "outputSource": "#main/s1/out2"
                }
            ],
            "requirements": [
                {
                    "types": [
                        {
                            "$import": "#recursive.yml/sample_meta"
                        },
                        {
                            "$import": "#recursive.yml/file_with_sample_meta"
                        },
                        {
                            "$import": "#recursive.yml/info_with_sample_meta"
                        },
                        {
                            "$import": "#array.yml/sample_meta2"
                        },
                        {
                            "$import": "#array.yml/study_meta"
                        },
                        {
                            "$import": "#array.yml/study_meta_too"
                        },
                        {
                            "$import": "#singletype.yml/simple_record"
                        },
                        {
                            "name": "#main/user_type1",
                            "type": "record",
                            "fields": [
                                {
                                    "name": "#main/user_type1/prop",
                                    "type": "string"
                                }
                            ]
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
