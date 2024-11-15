#!/usr/bin/env python

import argparse

import readfile
import rest_tool_functions

# dictionary for the allowed output formats
dict_output = {
    "cids": "xml",
    "aids": "xml",
    "sids": "xml",
    "description": "xml",
    "summary": "xml",
    "record": "csv",
    "classification": "xml",
    "targets": "txt",
    "xrefs": "txt",
    "synonyms": "txt",
    "property": "csv",
    "doseresponse": "csv",
}
check_for_id_type = ["cids", "aids", "sids"]

dic_key_value_type = {"assay": "AID", "compound": "CID", "substance": "SID"}
dic_key_value_operation = {"aids": "AID", "cids": "CID", "sids": "SID"}

post_id_types = ["inchi", "sdf", "smiles"]
id_dict = {"compound": "cid", "assay": "aid", "substance": "sid"}


def main(args):
    url = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/%s/%s/" % (args.type, args.id_type)
    # check if we are POST then skip this part otherwise insert the ids into the url
    if args.id_type not in post_id_types:
        if args.id_type == id_dict[args.type]:
            url += readfile.getListString(args) + "/"
        else:
            url += args.id_value + "/"

    url += args.operation + "/"
    if args.operation in ["target", "property", "xrefs"]:
        url += args.operation_value + "/"

    create_dict_tsv = False
    if args.operation == "xrefs":
        if "," in args.operation_value:
            url += "xml"
        else:
            url += "txt"
    else:
        if args.operation in check_for_id_type:
            # dont create dictionary if they are the same
            if dic_key_value_type[args.type] == dic_key_value_operation[args.operation]:
                url += "txt"
            else:
                url += "xml"
                create_dict_tsv = True
        else:
            url += dict_output[args.operation]
    if args.operation in check_for_id_type and args.id_type not in post_id_types:
        url += "?%s_type=%s" % (args.operation, args.ids_operation_type)
    print("The constructed REST URL is: %s" % url)

    if args.id_type in post_id_types:
        postfile = open(args.id_value, "r")
        post_value = postfile.read()
        post_dict = {args.id_type: post_value}
        # print(post_dict)
        readfile.store_result_post(url, post_dict, args.outfile)
    # check if we have to create a tsv file
    elif create_dict_tsv:
        key = dic_key_value_type[args.type]
        value = dic_key_value_operation[args.operation]
        dic = rest_tool_functions.get_dict_key_value(url, key, value)
        rest_tool_functions.write_to_sf(dic, args.outfile, "\t")
    else:
        readfile.store_result_get(url, args.outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--type", type=str, required=True, help="That you want BioAssay Compund ..."
    )
    parser.add_argument("--id-type", type=str, help="Specify the ID type")
    parser.add_argument(
        "--operation", type=str, required=True, help="Specify the operation"
    )
    parser.add_argument(
        "--operation-value",
        dest="operation_value",
        type=str,
        required=False,
        help="Specify the additional operation value",
    )
    parser.add_argument(
        "--xref-operation-value",
        dest="xref_operation_value",
        type=str,
        required=False,
        help="Specify the xref operation ",
    )
    parser.add_argument(
        "--ids-operation-type",
        dest="ids_operation_type",
        type=str,
        required=False,
        help="all inactive ...",
    )
    parser.add_argument(
        "--xref", dest="xref", type=str, help="use xref to identify the searched thing"
    )
    parser.add_argument(
        "--xref-value", dest="xref_value", type=str, help="Specify the xref"
    )
    parser.add_argument(
        "--property-value",
        dest="property_value",
        type=str,
        help="Specify the property value",
    )
    parser.add_argument(
        "--id-type-ff", dest="id_type_ff", type=str, help="file or field"
    )
    parser.add_argument(
        "--id-value", dest="id_value", type=str, required=True, help="Specify the id"
    )
    parser.add_argument(
        "--outfile",
        type=argparse.FileType("w"),
        required=True,
        help="Specify the output file",
    )

    args = parser.parse_args()
    main(args)
