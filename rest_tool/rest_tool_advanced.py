#!/usr/bin/env python

import sys, os
import argparse

import readfile

#dicitionary for the output format

dict_output={"cids" :"txt", "aids" : "txt", "sids" : "txt", "description": "xml", "summary" : "xml", "record" : "csv", "classification": "xml", "targets" : "txt", "xrefs" : "txt", "synonyms" : "txt", "property": "csv" }

#alles andere ist xml
check_for_id_type=["cids", "aids", "sids"]

id_dict={"compound": "cid", "assay": "aid", "substance" : "sid" }

def getListString(args):
    if args.id_type_ff == "file":
        #build comma list
        list_string=",".join(getListFromFile(open(args.id_value,"r")))
    else:
        list_string=args.id_value
    return list_string


def main(args):

    url="http://pubchem.ncbi.nlm.nih.gov/rest/pug/"+args.type+"/"
    url+=args.id_type+"/"
    if args.id_type ==id_dict[args.type]:
        url+=getListString(args)+"/"
    else:
        url+=args.id_value+"/"
    url+=args.operation+"/"
    if args.operation == "target" or args.operation == "property" or args.operation == "xrefs":
        url+=args.operation_value+"/"
    
    url+=dict_output[args.operation]
    if args.operation in check_for_id_type:
            url+="?"+args.operation+"_type="+args.ids_operation_type
    print(url)
    readfile.store_result(url, args.outfile)
    
    
    
    
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
