# A simple tool to connect to the AREsite server and retrieve feature
# information using the AREsite REST Interface.
# Parts of this code are from https://toolshed.g2.bx.psu.edu/repos/earlhaminst/ensembl_get_feature_info
import json
import optparse
import sys
import urllib.request, urllib.parse, urllib.error
import urllib.request, urllib.error, urllib.parse
import time
import requests
from six.moves.urllib.parse import urljoin

usage = "usage: %prog [options] arg1 arg2"
parser = optparse.OptionParser(usage=usage)
parser.add_option('-g', '--gene', help='Gene ID to search for')
parser.add_option('-m', '--motif', help='Motif to look for', default='ATTTA', type=str)
parser.add_option('-s', '--species', type='choice',
                  choices=['Homo_sapiens', 'Mus_musculus', 'Danio_rerio', 'Drosophila_melanogaster', 'Caenorhabditis_elegans'], default='Homo_sapiens',
                  help='Specify the species to investigate')
options, args = parser.parse_args()

if options.gene is None:
    raise Exception('- Specify the gene you want to look for!')

if "," in options.motif :
    raise Exception('- Please only search for single motifs at once')

class AREsiteRestClient(object):
    def __init__(self, server='http://rna.tbi.univie.ac.at/AREsite2/api/', reqs_per_sec=1):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, hdrs=None, params=None):
        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        if params:
            endpoint += '?' + urllib.parse.urlencode(params)

        data = None

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        try:
            request = urllib.request.Request(self.server + endpoint, headers=hdrs)
            response = urllib.request.urlopen(request)
            content = response.read().decode('utf-8')
            if content:
                data = json.loads(content)
            self.req_count += 1

        except urllib2.HTTPError as e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint, hdrs, params)
            else:
                sys.stderr.write('Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))

        return data

    def get_motifs(self, species, gene, motifs):
        query = str('?query={0}&species={1}&list={2}'.format(gene, species, motifs))
        if query:
            aresite = self.perform_rest_action(
                query
            )
            return aresite
        return None

def run(species, gene, motifs):
    client = AREsiteRestClient()
    aresite = client.get_motifs(species, gene, motifs)
    if aresite:

        mots        = aresite["exact_motifs"]
        starts      = aresite["motif_starts"]
        ends        = aresite["motif_ends"]
        chrs        = aresite["chromosomes"]
        strands     = aresite["strands"]
        transcripts = aresite["transcripts"]
        genes       = aresite["genes"]
        evh         = aresite["hur_evidence"]
        evt         = aresite["ttp_evidence"]
        eva         = aresite["auf_evidence"]
        anno        = aresite["annotation"]
        
        aresite = zip(chrs,starts,ends,mots,anno,strands,genes,transcripts,evh,evt,eva)

        def getKey(item):
            return item[1]

        aresite = sorted(aresite, key=getKey)
        
        for i in range(len(aresite)):
            print ("\t".join(aresite[i])+"\n")
                
            
if __name__ == '__main__':
    run(options.species, options.gene, options.motif)
