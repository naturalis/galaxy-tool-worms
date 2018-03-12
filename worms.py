#!/usr/bin/python
from SOAPpy import WSDL
import re
import sys, os, argparse
import json
import time
import requests

from fuzzywuzzy import fuzz
from fuzzywuzzy import process

parser = argparse.ArgumentParser(description = "worms api and gene checker")
requiredArguments = parser.add_argument_group('required arguments')
requiredArguments.add_argument('-i', '--input', metavar='input file', dest='input', type=str,
		    help='File with species names, on every line a new name', required=True)
requiredArguments.add_argument('-o', '--output', metavar='output file', dest='out', type=str,
		    help='A file with synonyms and sequence availability', required=True)
requiredArguments.add_argument('-m', '--mitofish', metavar='mitofish species input', dest='mitofish', type=str,
		    help='File containing mitofish species', required=True)
requiredArguments.add_argument('-k', '--klasse', metavar='klasse file input', dest='klasse', type=str,
		    help='File containing klasse export', required=True)
requiredArguments.add_argument('-s', '--12s', metavar='12S species file input', dest='genbank12s', type=str,
		    help='File containing species export from 12s selection', required=True)
requiredArguments.add_argument('-l', '--log', metavar='log file', dest='log', type=str,
		    help='Output file with errors', required=True)
requiredArguments.add_argument('-w', '--worms', metavar='worms check', dest='worms', type=str,
		    help='Yes if you want to check the worms database', required=True)
requiredArguments.add_argument('-f', '--fuzzy_klasse', metavar='fuzzy klasse check', dest='fuzzy_klasse', type=str,
		    help='Yes if you want to fuzzy check the klasse database', required=True)
requiredArguments.add_argument('-kf', '--fasta_klasse', metavar='fasta with klasse seqs', dest='fasta_klasse', type=str,
		    help='Yes if you want to output an fasta file with sequences from klasse', required=True)
requiredArguments.add_argument('-kfo', '--fasta_klasse_out', metavar='fasta klasse output', dest='fasta_klasse_out', type=str,
		    help='klasse fasta output', required=False, nargs='?', default="")
requiredArguments.add_argument('-bf', '--fasta_bold', metavar='fasta with bold seqs', dest='fasta_bold', type=str,
		    help='Yes if you want to output an fasta file with sequences from bold', required=True)
requiredArguments.add_argument('-bfo', '--fasta_bold_out', metavar='fasta bold output', dest='fasta_bold_out', type=str,
		    help='bold fatsa output', required=False, nargs='?', default="")

args = parser.parse_args()

def outputheader():
    with open(args.out, "a") as worms_output:
        header = ["#Input name", "#worms_scientificname", "#worms_valid_name", "#worms_valid_AphiaID", "#worms_isMarine", "#worms_isBrackish",
                  "#worms_Terrestrial", "#worms_isFreshwater", "#worms_status", "#worms_bold_id", "#worms_bold_link", "#worms_synonymes",
                  "#mitofish","#mitofish_hit", "#klasse_WScan", "#klasse_TWN", "#klasse_KRW", "#klasse_NZ", "#klasse_AphiaID", "#klasse_sequence", "#klasse_hit", "#12s_species",
                  "#12s_taxonid","#bold_sequence", "#bold_hit", "#klasse_fuzzy","#klasse_fuzzy_score"]
        worms_output.write("\t".join(header)+"\n")

def get_Records(name, offset_number):
    rec = requests.get('http://www.marinespecies.org/rest/AphiaRecordsByName/'+name,params={'like': "true", "marine_only": "false", "offset": offset_number}, allow_redirects=True)
    if rec.text and rec.text[0] != "<":
        a = json.loads(rec.text)
    else:
        a = False
    return(a)

def get_matchRecords(name):
    rec = requests.get('http://www.marinespecies.org/rest/AphiaRecordsByMatchNames?scientificnames%5B%5D='+name.replace(" ", "+"),params={"marine_only": "false"}, allow_redirects=True)
    if rec.text and rec.text[0] != "<":
        a = json.loads(rec.text)[0]
    else:
        a = False
    return(a)

def get_external_id(id, database):
    rec = requests.get('http://www.marinespecies.org/rest/AphiaExternalIDByAphiaID/' + str(id),params={"type": database}, allow_redirects=True)
    if rec.text and rec.text[0] != "<":
        a = json.loads(rec.text)[0]
        return a
    else:
        return ""

def get_sysnonyms_id(id):
    rec = requests.get('http://www.marinespecies.org/rest/AphiaSynonymsByAphiaID/' + str(id), allow_redirects=True)
    synNames = []
    if rec.text and rec.text[0] != "<":
        a = json.loads(rec.text)
        for synObject in a:
            synonym = synObject["scientificname"]
            if synonym:
                synNames.append(synonym.encode('utf-8'))
        return synNames
    else:
        return []

def make_mitofish_dict():
    with open(args.mitofish, "r") as mitofish_species:
        fishdict = {}
        for a in mitofish_species:
            fishdict[a.strip()] = a.strip()
    return fishdict

def mitofish_check(synonymsList, mitofishDict):
    mito = ["no", ""]
    for x in synonymsList:
        if x in mitofishDict:
            mito = ["yes", str(x)]
            break
    return mito

def make_klasse_dict():
    with open(args.klasse, "r") as klasse:
        klassedict = {}
        for line in klasse:
            line = map(str.strip, line.split("\t"))
            if line[-1]:
                seq = "yes"
            else:
                seq = "no"
            line.append(seq)



            if line[0].lower() in klassedict:
                if klassedict[line[0].lower()][-1] == "no":
                    klassedict[line[0].lower()] = line
            else:
                klassedict[line[0].lower()] = line
    return klassedict

def klasse_check(synonymsList, klasseDict):
    klasse = "\t".join(["", "", "", "", "", "", ""])
    for x in synonymsList:
        if x.lower() in klasseDict:
            if args.fasta_klasse == "yes":
                if klasseDict[x.lower()][-2]:
                    with open(args.fasta_klasse_out, "a") as klasse_fasta:
                        klasse_fasta.write(">" + x.lower() + "\n" + klasseDict[x.lower()][-2] + "\n")
            klasse = "\t".join(klasseDict[x.lower()][1:6])+"\t"+klasseDict[x.lower()][-1]+"\t"+x
            break
    return klasse

def make_12s_dict():
    with open(args.genbank12s, "r") as species:
        fishdict = {}
        for a in species:
            a = a.split("\t")
            fishdict[a[0].strip().lower()] = map(str.strip, a)
    return fishdict

def genbank12s_check(synonymsList, genbank12sDict):
    genbank = ["", ""]
    for x in synonymsList:
        if x.lower() in genbank12sDict:
            genbank = genbank12sDict[x.lower()]
            break
    return genbank

def bold_check(synonymsList):
    bold_result = ["no", ""]
    for x in synonymsList:
        if args.fasta_bold == "yes":
            rec = requests.get("http://www.boldsystems.org/index.php/API_Public/sequence?taxon="+x)
            if rec.text:
                bold_result = ["yes", x]
                with open(args.fasta_bold_out, "a") as bold_fasta:
                    bold_fasta.write(rec.text + "\n")
                break

    return bold_result

def make_klasse_species_list():
    with open(args.klasse, "r") as klasse:
        klasselist = []
        for line in klasse:
            line = line.split("\t")
            klasselist.append(line[0].lower())
    return klasselist

def klasse_fuzzy_check(synonymsList, klasselist):
    klasse_fuzzy = ["", ""]
    if args.fuzzy_klasse == "yes":
        fuzzy_hit = process.extractOne(synonymsList[0], klasselist)
        fuzzratioscore = fuzz.token_sort_ratio(synonymsList[0], fuzzy_hit[0])
        if fuzzratioscore > 80:
            klasse_fuzzy = [fuzzy_hit[0], str(fuzzratioscore)]

    return klasse_fuzzy

def get_all_worms_records(taxon_name):
    start = 1
    max_capacity = 50
    records = []
    record_data = {}
    a = get_Records(str(taxon_name).replace(" ","_"), start)
    if a:
        for i in a:
            records.append(i)
        while len(records) == max_capacity:
            start = start + 50
            max_capacity = max_capacity + 50
            b = get_Records(str(taxon_name.replace(" ","_")), start)
            if b:
                for i in b:
                    records.append(i)
        return records
    else:
        n = get_matchRecords(str(taxon_name))
        if n:
            for x in n:
                records.append(x)
            return records
        else:
            return None

def worms_check(taxon_name):
    hits = []
    if args.worms == "yes":
        records = get_all_worms_records(taxon_name)

        if records:
            for record in records:
                boldId = get_external_id(str(record["valid_AphiaID"]), "bold")
                boldLink = ""
                if boldId:
                    boldLink = "http://www.boldsystems.org/index.php/Taxbrowser_Taxonpage?taxid=" + str(boldId)
                synonymsList = [taxon_name, str(record["scientificname"]), str(record["valid_name"]), taxon_name]
                synonymsList.extend(get_sysnonyms_id(str(record["valid_AphiaID"])))

                record["input_name"] = taxon_name
                record["boldId"] = boldId
                record["boldLink"] = boldLink
                record["synonymsList"] = synonymsList
                hits.append(record)
        else:
            hits.append({"input_name":taxon_name, "synonymsList":[taxon_name]})
    else:
        hits.append({"input_name": taxon_name, "synonymsList": [taxon_name]})

    return hits

def worms(x, mitofishDict, klasseDict, genbank12sDict, klasse_species_list, output):
    tries = 3
    if x.strip():
        if x.strip()[-1] == ".":
            x = x.strip()[:-1]

        for i in range(tries):
            try:
                for hit in worms_check(x.strip()):
                    if len(hit) == 2:
                        tab = 12 * "\t"
                        output.write(hit["input_name"] + tab)
                    else:
                        output.write(hit["input_name"] + "\t" + str(hit["scientificname"]) + "\t" + str(
                            hit["valid_name"]) + "\t" + str(hit["valid_AphiaID"]) + "\t" + str(
                            hit["isMarine"]) + "\t" + str(hit["isBrackish"]) + "\t" + str(
                            hit["isTerrestrial"]) + "\t" + str(hit["isFreshwater"]) + "\t" + str(
                            hit["status"]).encode('utf-8') + "\t" + hit["boldId"].encode('utf-8') + "\t" + hit[
                                         "boldLink"].encode('utf-8') + "\t" + ", ".join(hit["synonymsList"]) + "\t")
                    output.write("\t".join(mitofish_check(hit["synonymsList"], mitofishDict)) + "\t")
                    output.write(klasse_check(hit["synonymsList"], klasseDict) + "\t")
                    output.write("\t".join(genbank12s_check(hit["synonymsList"], genbank12sDict)) + "\t")
                    output.write("\t".join(bold_check(hit["synonymsList"])) + "\t")
                    output.write("\t".join(klasse_fuzzy_check(hit["synonymsList"], klasse_species_list)))
                    output.write("\n")
            except:
                if i < tries - 1:
                    continue
                else:
                    with open(args.log, "a") as error_log:
                        e = str(sys.exc_info()[0]) + " " + str(sys.exc_info()[1])
                        error_log.write("error: " + x.strip() + "\t" + str(e) + "\n")
            break


def check_input():
    mitofishDict = make_mitofish_dict()
    klasseDict = make_klasse_dict()
    genbank12sDict = make_12s_dict()
    klasse_species_list = make_klasse_species_list()

    with open(args.input, "r") as fish_list, open(args.out, "a") as output:
         for x in fish_list:
             worms(x, mitofishDict, klasseDict, genbank12sDict, klasse_species_list, output)


def main():
    outputheader()
    check_input()

if __name__ == '__main__':
    main()

