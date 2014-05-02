#!/usr/bin/env python

import sys, os
import argparse

import readfile

#dicitionary for the output format

dict_output={"cids" :"txt", "aids" : "txt", "sids" : "txt", "description": "xml", "summary" : "xml", "record" : "csv", "classification": "xml", "targets" : "txt", "xrefs" : "txt", "synonyms" : "txt", "property": "csv", "doseresponse" : "csv" }

#alles andere ist xml
check_for_id_type=["cids", "aids", "sids"]

post_id_types=["inchi", "sdf", "smiles"]

id_dict={"compound": "cid", "assay": "aid", "substance" : "sid" }

def getListString(args):
    if args.id_type_ff == "file":
        #build comma list
        list_string=",".join(getListFromFile(open(args.id_value,"r")))
    else:
        print (args.id_value)
        list_string=args.id_value.strip().replace("__cr____cn__", ",")
    return list_string

    
def main(args):
    url="http://pubchem.ncbi.nlm.nih.gov/rest/pug/"+args.type+"/"
    url+=args.id_type+"/"
    # check if we are post then skip this part otherwise put the ids in the url
    if not args.id_type in post_id_types:
        if args.id_type ==id_dict[args.type]:
            url+=getListString(args)+"/"
        else:
            url+=args.id_value+"/"

    url+=args.operation+"/"
    if args.operation == "target" or args.operation == "property" or args.operation == "xrefs":
        url+=args.operation_value+"/"
    
    if args.operation == "xrefs":
        if "," in args.operation_value:
            url+="xml"
        else:
            url+="txt"
    else:
        url+=dict_output[args.operation]
    if args.operation in check_for_id_type and args.id_type not in post_id_types:
            url+="?%s_type=%s" % (args.operation, args.ids_operation_type)
    print(url)
    if args.id_type in post_id_types:
        postfile=open(args.id_value,"r")
        post_value=postfile.read()
        post_dict={args.id_type : post_value}
        print(post_dict)
        readfile.store_result_post(url, post_dict, args.outfile)
    else:
        readfile.store_result_get(url, args.outfile)
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--type', type=str, required=True,
        help="That you want BioAssay Compund ...")
    parser.add_argument('--id-type', type=str,
        help="Specify the ID type")
    parser.add_argument('--operation', type=str, required=True,
        help="Specify the operation")
    parser.add_argument('--operation-value', dest="operation_value", type=str, required=False,
        help="Specify the additional operation value")
    parser.add_argument('--xref-operation-value', dest="xref_operation_value", type=str, required=False,
        help="Specify the xref operation ")
    parser.add_argument('--ids-operation-type', dest="ids_operation_type", type=str, required=False,
        help="all inactive ...")
    parser.add_argument('--xref', dest="xref", type=str,
        help="use xref to identify the searched thing")
    parser.add_argument('--xref-value', dest="xref_value", type=str,
        help="Specify the xref")
    parser.add_argument('--property-value', dest="property_value", type=str,
        help="Specify the property value")
    parser.add_argument('--id-type-ff', dest="id_type_ff", type=str,
        help="file or field")
    parser.add_argument('--id-value', dest="id_value", type=str, required=True,
        help="Specify the id")
    parser.add_argument('--outfile', type=argparse.FileType('w'), required=True,
        help="Specify the output file")


    if len(sys.argv) < 8:
        print "Too few arguments..."
        parser.print_help()
        exit(1)
    args = parser.parse_args()
    main( args )
