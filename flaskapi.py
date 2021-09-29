#!/usr/bin/python3.6
# Flask api for regsp
# Author: Girum Fitihamlak Ejigu
# Last modified on September 2021

from flask import Flask, request, redirect, jsonify, send_from_directory
from flask_cors import CORS
from flask.globals import current_app
from werkzeug.utils import secure_filename
import os, json
from Bio import SeqIO
from Bio import Entrez
import sqlite3
import glob
from time import strftime
import time
import re
import shutil,subprocess
import threading
from gsp import READ_SEQ, FIND_SUBREGION_CTG, HANDLE_RAW_DATA
# File upload folder
Entrez.email = "email@domain.com"
Entrez.api_key= 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
UPLOAD_FOLDER = '/home/user/regsp/Files'
ALLOWED_EXTENSIONS = {"txt", "fasta", 'fa'}

def allowed_files(filename):
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXTENSIONS


app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
CORS(app)

tasklists={}
@app.route('/api', methods=["GET","POST"])
def sendPostt():
    if request.method == "POST":
        #print (request.form, flush=True)
        ipaddr=request.remote_addr
        #print(os.getcwd())
        if not os.path.exists("/home/user/regsp/Files/"+ipaddr):
            os.mkdir(ipaddr)
        accesion = request.form.get("accession")
        max_allowed_hits = int(request.form.get("maxb"))
        minDBcount = int(request.form.get("ming"))
        query = request.files["query"]
        selected = json.loads(request.form.get("selectedAcc"))
        kmersize=request.form.get("kmer")
        if(request.form.get("db")=='null'):
            userref=None
            userdb=None
        else:
            userref=request.files["ref"]
            userdb=request.files["db"]
        #print(userdb)
        
        if query and allowed_files(query.filename):
            filename = secure_filename(query.filename)
            query.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            if(userdb is not None):
                userref.save(os.path.join(app.config['UPLOAD_FOLDER'], userref.filename))
                userdb.save(os.path.join(app.config['UPLOAD_FOLDER'], userdb.filename))

            # take accecceion number from text box or from list search
            try:
                if(userdb is not None):
                    userreflist=[]
                    #userref.save(os.path.join(app.config['UPLOAD_FOLDER'], userref.filename))
                    with open(userref.filename, "r") as fa:
                        lines=fa.read().split('>')
                        lines = lines[1:]
                        lines=['>'+ seq for seq in lines]
                        for name in lines:
                            #Extracting sequence Id to use it for file name
                            file_nameR=name.split('\n')[0][1:]  
                            out_file=open(file_nameR+".fasta", "w")
                            out_file.write(name)
                            out_file.close()
                            userreflist.append(os.path.abspath(file_nameR+".fasta"))
                    refFileName=userreflist
                    #print(refFileName)

                    userdblist=[]
                    #userdb.save(os.path.join(app.config['UPLOAD_FOLDER'], userdb.filename))
                    with open(userdb.filename, "r") as f:
                        buff = []
                        i = 1
                        file_name=userdb.filename.split(".")[0].split()[0]
                        for line in f:
                            if line.strip():
                                buff.append(line)
                            if line.strip()[0:3] == "***":
                                output = open(file_name+'%d.fasta' % i,'w')
                                output.write(''.join(buff))
                                output.close()
                                userdblist.append(os.path.abspath(file_name+"%d.fasta" %i))
                                i+=1
                                buff = []
                        output = open(file_name+'%d.fasta' % i,'w')
                        output.write(''.join(buff))
                        output.close()
                        userdblist.append(os.path.abspath(file_name+"%d.fasta" %i))
                        dbFileName,emtyCDS = checkCDS(userdblist)
                        if len(dbFileName)==0:
                            return jsonify( data = "No_CDS")
                        #print(dbFileName)
                    

                elif (not accesion or accesion=="null"):
                    #print (selected[0])
                    codingsequences  = downloadCodingsequences(selected)
                    dbFileName, emtyCDS =checkCDS(codingsequences)
                    if len(dbFileName)==0:
                        return jsonify( data = "No_CDS")
                    refFileName = downloadrefsequences(selected)

                else:
                    acclist=accesion.split(",") # get list of accesions without space
                    #print(acclist)
                    codingsequences  = downloadCodingsequences(acclist)
                    dbFileName,emtyCDS = checkCDS(codingsequences)
                    if len(dbFileName)==0:
                        return jsonify( data = "No_CDS")
                    refFileName = downloadrefsequences(acclist)
            except:
                return jsonify( data = "filerror")
            appl=current_app._get_current_object()
            th=threading.Thread(target=reGSP, args=(appl,filename,dbFileName,refFileName,max_allowed_hits,minDBcount,kmersize,ipaddr))
            th.start()
            if ipaddr not in tasklists.keys():
                tasklists[ipaddr]=[]
            tasklists[ipaddr].append(th)
            return jsonify( data = "Submitted")
        else:
            return jsonify( data = "filerror")


def reGSP(app,filename,dbFileName, refFileName,max_allowed_hits,minDBcount,kmersize,ipaddr):
    #print ("File Uploaded")
    with app.app_context():   
        seqFileName = filename
        seqFileNamecPlot=os.path.abspath(seqFileName)
        refFileName2cPlot=joinfiles(refFileName)
        #print(seqFileNamecPlot,refFileName2cPlot)
        timest=time.strftime("%Y-%m-%dT%H:%M:%S")
        outputFileName = filename.split(".")[0]+"@"+timest

        #PRINT_STR("=> Reading a query sequence...\n")
        seq_db, seq_db_idx = READ_SEQ(seqFileName)

        #PRINT_STR("=> Calculating...\n")
        
        dirr="/home/user/regsp/Files/"+ipaddr
        
        raw_hits, raw_data = FIND_SUBREGION_CTG(seq_db, seq_db_idx, dbFileName, outputFileName,max_allowed_hits)
        
        HANDLE_RAW_DATA(seq_db, seq_db_idx, raw_hits, raw_data, refFileName, outputFileName, minDBcount,dirr)
        cplot='/home/user/regsp/cPlot.py'
        cmd="python %s %s %s %s %s %s" % (cplot,refFileName2cPlot,seqFileNamecPlot,dirr,outputFileName,kmersize)
        cmd_out=subprocess.getoutput(cmd)
        if(whitespace_only(os.path.join(dirr,outputFileName+".plot.seq"))):
            return jsonify(data = "No_Plot"+"-"+outputFileName)

        return jsonify( data = outputFileName)

@app.route('/inprog', methods =["GET"])
def inprogress():
    ipad=request.remote_addr
    tasklist={}
    #print('tasklists----->', len(tasklists))
    if ipad in tasklists.keys():
        tasklist[ipad] = [t for t in tasklists[ipad] if t.is_alive()]
        #print("tasklist IP ---->",len(tasklists[ipad]))
    else:
        tasklist[ipad]=[]
    #print("Getting update")
    if tasklist[ipad]==[]:
        if ipad in tasklists.keys():
            del(tasklists[ipad])
        return jsonify(data=0)
    #print("length *** ",len(tasklist[ipad]))
    return jsonify(data=len(tasklist[ipad]))

@app.route('/uploads/<path:filename>')
def download_file(filename):
    ipaddress=request.remote_addr
    directory="/home/user/regsp/Files/"+ipaddress
    return send_from_directory(directory,
                               filename, as_attachment=True)

def downloadrefsequences(accessionno):
    filenames=[]
    for idx in accessionno:
        filename= idx+".fasta"
        if not os.path.isfile(filename):
            net_handle = Entrez.efetch(
                db="nucleotide", id=idx, rettype="fasta", retmode="text"
            )
            out_handle = open(filename, "w")
            out_handle.write(net_handle.read())
            out_handle.close()
            net_handle.close()
        filenames.append(os.path.abspath(filename))
    return filenames

def downloadCodingsequences(accessionno):
    filenames=[]
    for idx in accessionno:
        filename = idx+"CDS.fasta"
        if not os.path.isfile(filename):		
            net_handle = Entrez.efetch(
                db="nucleotide", id=idx, rettype="fasta_cds_aa", retmode="text"
            )
            #if(net_handle.read()==""):
            #    return "No_CDS"
            out_handle = open(filename, "w")
            out_handle.write(net_handle.read())
            out_handle.close()
            net_handle.close()
        filenames.append(os.path.abspath(filename))
    return filenames

def whitespace_only(file):
    content = open(file, 'r').read()
    if re.search(r'^\s*$', content):
        return True
def checkCDS(cdsList):
    coding=[]
    ncoding=[]
    for idx in cdsList:
        if whitespace_only(idx):
            ncoding.append(idx)
        else:
            coding.append(idx)
    return coding, ncoding

def joinfiles(filelist):
    with open('Ref_file.fasta','wb') as wfd:
        for f in filelist:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)
    return os.path.abspath(wfd.name)

@app.route('/checkAccession', methods =["GET","POST"])
def showSearchedReferenceFile():
    sterm=request.form.get("searchTerm")
    accessionSearch=request.form.get("accSearch")
    if (accessionSearch=="true"):
        invalidAccs, validAccs = validAccessionSearch(sterm)
        handle= Entrez.efetch(db="nucleotide", id=validAccs, rettype="gb", retmode="text")
        result=handle.read()
        handle.close()
        gbfiles=result.split("//\n")  ## split genebank text files to a list
        data=[{"accession":a,"file":b} for a,b in zip(validAccs,gbfiles)] ## create a dictonary list with accesion and genebank file
        return jsonify(data = data, notfound=invalidAccs)
    else:
        accession=searchNCBI(sterm)
        if (accession==[]):
            return jsonify(data = None)
        else:
            handle= Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            result=handle.read()
            handle.close()
            gbfiles=result.split("//\n")  ## split genebank text files to a list
            data=[{"accession":a,"file":b} for a,b in zip(accession,gbfiles)] ## create a dictonary list with accesion and genebank file
            return jsonify(data = data)


def searchNCBI(sterm):
    searchTerm ="("+sterm+"[ACCN] AND srcdb_refseq[prop]) OR ("+sterm+"[TI] AND srcdb_refseq[prop])" 
    handle = Entrez.esearch(db="nucleotide", retmax=20, term=searchTerm, idtype="acc")
    result=Entrez.read(handle)
    handle.close
    #print(result["Count"])
    return result["IdList"]

def validAccessionSearch(sterm):
    invalidAcc=[]
    validAcc=[]
    accessionList=sterm.split(",")
    for x in accessionList:
        if (searchNCBI(x)==[]):
            invalidAcc.append(x)
        else:
            validAcc.append(x)
    return invalidAcc, validAcc

@app.route('/searchTree/<searchdbTerm>')
def searchDatabase(searchdbTerm):
    conn = sqlite3.connect('Organelles.db')
    curr = conn.cursor()
    query="""select json_object('name', A3.supergroup, 'children', json_group_array(json(A3.json_obj3))) from
            (select A2.supergroup, json_object('name', A2.kingdom, 'children', json_group_array(json(A2.json_obj2)))as json_obj3 from 
                (select A1.supergroup, A1.kingdom, json_object('name', A1.family, 'children', json_group_array(json(A1.json_obj1))) as json_obj2 from 
                    (select supergroup, kingdom, family, json_object('name', genus, 'children', json_group_array(json_object('name',organismName, 'accession',accession, 'type',type)))as json_obj1
                    from organelle where organismName like'%{}%'
                    group by supergroup, kingdom, family, genus) as A1
                group by A1.supergroup, A1.kingdom, A1.family) as A2
            group by A2.supergroup, A2.kingdom) as A3
        group by A3.supergroup """.format(searchdbTerm)
    curr.execute(query)
    result=curr.fetchall()
    if(result==[]):
        conn.commit()
        conn.close()
        return jsonify(data=None)
    else:
        searchdbResult="[{}]".format(result[0][0])
        conn.commit()
        conn.close()

        return jsonify(data=searchdbResult)

@app.route('/genbank/<accession>')
def getGenebank(accession):
    net_handle = Entrez.efetch(
            db="nucleotide", id=accession, rettype="gb", retmode="text"
        ) 
    gbfile=net_handle.read()
    net_handle.close()
    return jsonify(data = gbfile)

@app.route('/addHigher', methods =["GET","POST"])
def addhiger():
    query=request.form.get("query")
    conn = sqlite3.connect('Organelles.db')
    curr = conn.cursor()
    curr.execute(query)
    result=curr.fetchall()
    keys=['name','accession','type']
    d=[dict(zip(keys,l))for l in result]
    addHigerNodeResult=json.dumps(d)
    conn.commit()
    conn.close()

    return jsonify(data=addHigerNodeResult)

@app.route('/history', methods =["GET","POST"])
def checkHistory():
    ipaddress=request.remote_addr
    if not os.path.exists("/home/user/regsp/Files/"+ipaddress):
        return jsonify(data=None)
    else:
        filepath="/home/user/regsp/Files/"+ipaddress+"/"
        seqf=glob.glob(filepath+"*.seq")
        seqf.sort(key=os.path.getmtime)
        lastseq=[os.path.basename(i) for i in seqf[::-1]]
        sortedf=glob.glob(filepath+"*.sorted")
        sortedf.sort(key=os.path.getmtime)
        lastsorted=[os.path.basename(i) for i in sortedf[::-1]]
        regspf=glob.glob(filepath+"*.regsp")
        regspf.sort(key=os.path.getmtime)
        lastregsp=[os.path.basename(i) for i in regspf[::-1]]
        plotf=glob.glob(filepath+"*.plot.pdf")
        plotf.sort(key=os.path.getmtime)
        lastplot=[os.path.basename(i) for i in plotf[::-1]]
        cPlotf=glob.glob(filepath+"*.cPlot.pdf")
        cPlotf.sort(key=os.path.getmtime)
        lastcPlot=[os.path.basename(i) for i in cPlotf[::-1]]

        datetime=[i.split("@")[1].split(".")[0] for i in lastsorted]
        history=[{"date":a, "query":b,"sorted":c,"plot":d,"cPlot":e, "regsp":f} for a, b, c, d, e, f in zip(datetime, lastseq, lastsorted, lastplot,lastcPlot,lastregsp)]
        return jsonify(data=history)



if __name__ == "__main__":
    app.run()
