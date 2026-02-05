#!/usr/bin/env python3
# Vendored from MAD v2.2 (27-Mar-2018) by Tria, Landan & Dagan.
# Source: https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen
# This file is included unmodified except for this header comment.
#=====================
"""
MAD phylogenetic rooting
       - Please cite DOI:10.1038/s41559-017-0193

Usage:  
        mad filename [-mnsptfgvh]
    
    Where the file contains tree(s) in NEWICK format.
    Rooted tree(s) will be written to 'filename.rooted',
    rooting statistics to screen.

    Flags:
       -m: Multiple trees in input file, output file is verbose and
           contains '>' (info) and '<' (error) lines in addition to 
           NEWICK strings. 
           Default allows exactly one tree, and output is pure NEWICK.
       -n: Like -m, but pure newick. Only one rooted tree per input tree,
           errors as a lone ';' in a line.
       -s: Statistics - append MAD statistics to output NEWICK strings.
       -p: Polytomies - flat polytomies in rooted NEWICK (if present).         
           Default is a binary tree with polytomies arbitrarily resolved
           into zero-length branches.       
       -t: Tiny - retain branch lengths smaller than 10^-6.                
           Default contracts these to 0.0, thereby creating polytomies.  
  -f | -g: Figtree - rooted trees formatted for viewing AD scores in figtree.
           Important: after loading in figtree, please check the option 
           'Appearance>Gradient' manually, otherwise the branch colors will
           be misleading. -f reports ancestor deviations (AD) only for nodes,
           while -g reports also within-branch maximizing positions and AD values.           
       -v: Version.       
       -h: Help.
        
Please report bugs to giddy.landan@gmail.com .
"""
#v2.2 27-Mar-2018 
#=====================
#---- init
from sys import argv, stderr, version, exit
from re import sub
fstr="".join([a.upper() for a in argv if a[0]=='-'])
flags={f: fstr.count(f) for f in "MNSPTFGVHD"}
if flags['V']: 
    v=sub('\s+[^\[]+',' ',version,count=1)
    print('mad 2.2\npython '+v)
    exit()
if flags['G'] and flags['F']:
    print('Flags -g and -f are mutually exclusive.')
    exit()
if flags['M'] and flags['N']:
    print('Flags -m and -n are mutually exclusive.')
    exit()
if len(argv)<2 or flags['H']:
    print(__doc__)
    exit()
fstr=[a for a in argv if a[0]!='-']
if len(fstr)!=2:
    exit('Expecting exactly one filename')
print("\nMAD phylogenetic rooting")
fn=fstr[1]
#-------- consts
if flags['T']: minlen=1e-15
else:          minlen=1e-6
madtol=1.0001
gdbug_level=flags['D']
at="#^$%"
if flags['M'] or flags['N']:  err="Error: "
else:           err="Error analyzing file '"+fn+"':\n\t"
bug=err+"Oops... Please report this error to giddy.landan@gmail.com ."
class madError(Exception):
    pass
#------- short circut init - deffered imports
try:
    with open(fn) as x: nwkstr=x.read()
except FileNotFoundError:           exit("\nFile not found: "+fn)
except:                             exit("\nError reading file: "+fn)
blackspace=str.maketrans('', '', '\t\n\r\f\v')
mnwk=nwkstr.translate(blackspace).split(';')
if mnwk[-1]=='': mnwk.pop()
nnwk=len(mnwk)
if nnwk==0:                    exit(err+
                                    "No NEWICK string found.")
if nnwk>1 and not flags['M'] and not flags['N']:  exit(err+
                                    "Expecting exactly one string terminated with ';'. (For multiple inputs, use -m.)")
#---
ofn=fn+".rooted"
try:
    fh=open(ofn,'w')
except:
    exit("\nError opening output file: "+ofn)
#---------------- 
#------ imports
import re
import numpy as np
from copy import deepcopy
from statistics import stdev, mean
from inspect import stack
#----- functions
#========================================
#===== write to output file
def writeout(s):
    global fh,ofn
    try:    
        fh.write(s+"\n")
    except:
        fh.close() 
        exit("\nError writing to file: "+ofn)
#<-- writeout
#========================================

#========================================
#===== deal outputs
def madlog(s):
    global fh,ofn,gdbug_level,flags
    if s[0]=='>':           pref=''
    elif s[:5]=='Error':    pref='<<< '
    else:                   pref='>> '
    if gdbug_level<1:       
        pos=''
    else:
        sat=stack()[1][1:3]
        pos="{}:{} ".format(sat[0],sat[1])
    if pref=='' and flags['M']:
        writeout(s.strip())
    else:
        s1=re.sub('(?m)^',pref+pos,s.strip())
        if s1[0]=='<': s1+="\n"    
        if flags['M']: writeout(s1)
        print(s1)
    if s[:5]=='Error':
        if flags['N']:
            writeout(';')
        raise madError
#<-- madlog
#========================================

#========================================
#===== debug messages
def gdbug(lev,*args):
    global gdbug_level
    if lev>gdbug_level:
        return
    sat=stack()[1][1:3]
    print("---\n",file=stderr)
    [print(x,file=stderr) for x in args]
    print("<-- {}:{}".format(sat[0],sat[1]),file=stderr)   
#<-- gdbug
#========================================

#========================================
#===== branch lengths conversion and sanity
def str2len(s):
    try:
        x=float(s)
    except ValueError:      madlog(err+"Corrupt NEWICK format: invalid branch length (...'"+s+"'...).")
    if not np.isfinite(x):  madlog(err+"Cowardly refusing to root trees with infinite or undefined branch lengths (...'"+s+"'...).") 
    return x
#<-- str2len
#========================================

#========================================
#======= convert NEWICK string, fills: labs,blen,nodes,trip + various scalars
def nwk2tree(nwkstr):
    nc=[ nwkstr.count(x) for x in ';,:()']
    if nc[0]!=1 or nwkstr[-1]!=';':     madlog(err+
                                        "Corrupt NEWICK format - expecting exactly one string terminated with ';'.")
    if nc[3]!=nc[4]:                    madlog(err+
                                        "Corrupt NEWICK format - unbalanced ().")
    notu=nc[1]+1
    if notu<3:                          madlog(err+
                                        "Cowardly refusing to root trees with less than 3 OTUs.") 
    if nc[2]>2*notu-2:                  madlog(err+
                                        "Corrupt NEWICK format - too many ':' | too few ','.") 
    if nc[3]>notu-1:                    madlog(err+
                                        "Corrupt NEWICK format - too many '()' | too few ','.") 
    nnode=notu*2-1
    rootnode=nnode-1
    nodes=np.full((nnode,3),-1,np.int16)
    trip=np.full((nnode,nnode),2,np.int8)
    blen=[0]*nnode
    labs=['']*nnode
    cld_splt=re.compile('\(([^()]*)\)')
    node_splt=re.compile('[:@]')
    iotu,tlen,ntiny,prec=[0]*4
    badbsp=[]
    inode=notu
    s2=nwkstr.replace("@",at).strip()
    while s2[0]=="(" :
        c3=cld_splt.split(s2,1)
        if len(c3)!=3:                  madlog(err+
                                        "Corrupt NEWICK format: () do not balance.") 
        b2=c3[1].split(',')
        if len(b2)<2:                   madlog(err+
                                        "Corrupt NEWICK format: empty, singleton or unbalanced group ().") 
        elif len(b2)>2:
            b2[1]="({},{}):0.0".format(b2[0],b2[1])
            z=str.join(',',b2[1:])
            s2="{}({}){}".format(c3[0],z,c3[2])        
            continue
        for i in [0,1]:
            b=b2[i]
            if b[0]!="@":
                z=b.split(":")
                if len(z)!=2:           madlog(err+
                                        "Corrupt NEWICK format: malformed 'node_label:branch_length' clouse (...'"+b+"'...).") 
                if labs.count(z[0])>0:  madlog(err+
                                        "Cowardly refusing to root trees with duplicate OTU names ('"+z[0]+"').") 
                b="@{}@{}".format(iotu,b)
                trip[iotu,iotu]=-1
                iotu+=1
            ll=node_splt.split(b)
            if len(ll)!=4:              madlog(err+
                                        "Corrupt NEWICK format: malformed 'node_label:branch_length' clouse (...'"+ll[2].replace(at,"@")+"'...).") 
            jnode=int(ll[1])
            labs[jnode]=ll[2]
            x=str2len(ll[3])
            if abs(x)<minlen and x!=0 and not flags['T']:
                #gdbug(3,ll[3],x)
                x=0
                ntiny+=1
            if x<0:                     madlog(err+
                                        "Cowardly refusing to root trees with negative branch lengths.") 
            blen[jnode]=x
            tlen+=x
            prec=max(prec,len(ll[3]))
            nodes[jnode,2]=inode
            nodes[inode,i]=jnode
            trip[inode,trip[jnode,:]<2]=i
        #<--i 0:1      
        trip[inode,inode]=-1
        s2="{}@{}@{}".format(c3[0],inode,c3[2])        
        inode+=1
    #<--inode loop            
    n=sum([c3[2].count(x) for x in '(),:'])
    if n>0:                             madlog(err+
                                        "Corrupt NEWICK format: () balanced out before end of string, tail contains more '(),:'.\n"+
                                        "\t(Possible cause - missing ';' between two trees).") 
    #gdbug(1,[notu,iotu,nnode,inode,tlen,prec,ntiny],badbsp)
    if inode!=nnode or iotu!=notu:      
        #gdbug(2,s2,nwkstr,nodes)
        print(s2)
        exit()
        madlog(bug+' (Unidentified parsing error.)\nInput string: '+nwkstr)
    if  ntiny>0 and not flags['T']: madlog("Warning: {} tiny branch lengths (<10^-6) were converted to 0. (Override with '-t'.)".format(ntiny))
    n=sum([x>0 for x in blen])
    if n<3:                             madlog(err+
                                        "Cowardly refusing to root trees with less than 3 positive branch lengths.") 
    if tlen==0:                         madlog(err+
                                        "Cowardly refusing to root zero-lengthed trees.")    
    slen=np.sort(blen)
    n=sum(np.diff(slen[slen>0])==0)
    if n>0: madlog("Warning: Trees with repeating branch lengths are suspicious ({} repeating values).".format(n))
    if len(badbsp)>0: madlog("Warning: Non-numeric or negative bootstrap values ignored. ({}).".format(badbsp))
    #gdbug(2,nodes,blen,trip,labs)
    return(notu,nnode,rootnode,nodes,blen,trip,labs,prec,tlen)
#<--nwk2tree()
#========================================

            
#========================================
#-------- Pre-processing:
def pwdist_preproc():
    #======= pairwise distances
    dist=np.full((nnode,nnode),0.0,np.float64)
    for i in range(notu,nnode):
        for j in [0,1]:
            b=nodes[i,j]
            for k in np.flatnonzero(trip[b,:]<2):
                dist[i,k]=dist[b,k]+blen[b]
                dist[k,i]=dist[i,k]
        for k0 in np.flatnonzero(trip[i,:]==0):
            for k1 in np.flatnonzero(trip[i,:]==1):
                dist[k0,k1]=dist[i,k0]+dist[i,k1]
                dist[k1,k0]=dist[k0,k1]
    #======= tip polytomies
    etrip=deepcopy(trip)
    n2n=list(range(nnode))
    ntipp=0            
    for i in range(nnode):
        if n2n[i]<i:
            continue
        jj=np.flatnonzero(dist[i,:]==0)
        if len(jj)>1:
            for j in jj[1:]:
                n2n[j]=i
                if i<notu and j<notu:
                    etrip[j,:]=69
                    etrip[:,j]=69
                    ntipp+=1   
    enotu=notu-ntipp            
    npair=enotu*(enotu-1)/2          
    #np.set_printoptions(threshold=np.inf)                        
    #gdbug(1,notu,ntipp,enotu,npair,n2n)            
    if ntipp>0: madlog("Warning: Squeezing tip polytomies ({} OTUs, {} redundant tips, {} effective OTUs).".format(notu,ntipp,enotu))
    if enotu<3: madlog(err+"Cowardly refusing to root trees with less than 3 effctive OTUs.") 
    #======= preproc in lists
    kij=[[]]*4
    dij=[[]]*nnode
    for i in range(nnode):
        for j in [0,1,2]:
            kij[j]=np.flatnonzero(etrip[i,:notu]==j)
        kij[3]=np.flatnonzero(etrip[i,:notu]<2)
        dij[i]=list([[]]*4)
        for j in [0,1,2,3]:
            dij[i][j]=dist[i,kij[j]].flatten().tolist()
    #gdbug(2,dij,"\n",npair)            
    return(enotu,npair,dist,dij,n2n)    
#<-- pwdist_preproc()
#========================================

#========================================
#-------- branch AD
def ancestor_deviations():
    global dist
    #-------- node deviations triplets
    dsum=[[]]*nnode
    for i in range(nnode):
        dsum[i]=[0]*3
        if i<notu: continue
        for j in [0,1,2]:
            k1=(j+1)%3
            k2=(j+2)%3
            ndev=0
            for dik1 in dij[i][k1]:
                for dik2 in dij[i][k2]:
                    d=dik1+dik2
                    if d>0:
                        ndev+=(2*dik1/d - 1)**2
            dsum[i][j]=ndev
    #<-- for i
    #-------- node ad
    nad=[-1]*nnode
    #nccv=[-1]*nnode
    #ndepth=[-1]*nnode
    #r2t=[]
    ibrn=range(notu,nnode)
    for i in range(nnode):
        if n2n[i]<i:
            nad[i]=nad[n2n[i]]
            continue
        if i<notu:
            nomin=enotu-1 
        else:
            nomin=sum(dsum[i])
        for k in ibrn: 
            if k!=i:
                nomin+=dsum[k][trip[k,i]]
        #gdbug(2,nomin,npair) 
        nad[i]=(nomin/npair)**0.5
        if nad[i]>1.0:
            madlog(bug+' (Node AD out of range.)\nInput tree is: '+nwkstr)           
        #r2t=dij[i][3]+dij[i][2]
        #nccv[i]=stdev(r2t)/mean(r2t)
        #ndepth[i]=max(r2t)        
    #-------- transversing pairs and branch ad
    ad=[10**6]*nnode
    rlen=[-1]*nnode
    rlen2=[-1]*nnode
    ccv=[-1]*nnode
    depth=[-1]*nnode
    polyroots=[-1]*nnode
    #gdbug(1,blen)
    for i in range(nnode-1):
        if blen[i]==0: 
            rlen[i]=0
            rlen2[i]=0
            ad[i]=nad[i]+1
            continue
        j=nodes[i,2]
        jj=trip[j,i]
        denom=0
        nomin=0    
        for dik1 in dij[i][3]:
            for dik2 in dij[i][2]:
                d=dik1+dik2 
                d2=d**-2
                denom+=d2
                nomin+=d2*(dik2-dik1)
        rho=nomin/(2*denom)
        r1=min(max(0,rho),blen[i])
        r2=blen[i]-r1
        rn=-1
        if r1<minlen and r1<r2:
            r1=0
            r2=blen[i]
            rn=i
        elif r2<minlen and r2<r1:
            r1=blen[i]
            r2=0
            rn=j
        rlen[i]=r1
        rlen2[i]=r2
        if rn>-1 and polyroots[n2n[rn]]>-1:
            ad[i]=nad[rn]+1            
            continue
        nomin=0    
        for dik1 in dij[i][3]:
            for dik2 in dij[i][2]:            
                d=dik1+dik2 
                nomin+=(2*(dik1+r1)/d - 1)**2
        for k in ibrn: 
            if trip[i,k]==2:
                nomin+=dsum[k][trip[k,i]]
            else:
                nomin+=dsum[k][trip[k,j]]
        ad[i]=(nomin/npair)**0.5
        if ad[i]>1.0:
            madlog(bug+' (AD out of range.)\nInput tree is: '+nwkstr)           
        r2t=[ d+r1 for d in dij[i][3]]+[ d-r1 for d in dij[i][2]]
        ccv[i]=stdev(r2t)/mean(r2t)
        depth[i]=max(r2t)
        if rn>-1:
            polyroots[n2n[rn]]=i
    #gdbug(1,ad,nad,rlen,rlen2,ccv,polyroots)
    return(ad,nad,rlen,rlen2,ccv,depth)
#<-- ancestor_deviations()            
#========================================

#========================================
#======= convert to NEWICK string
def tree2nwk():
    global labs
    nstck=[rootnode,nodes[rootnode,0],nodes[rootnode,2]]
    nwk=[[]]*nnode
    nwk[nstck[1]]=[rootnode]
    nwk[nstck[2]]=[rootnode]
    #gdbug(2,"yy {} {} {} yy".format(rootnode,nwk,nstck))
    blens=[prec.format(x) for x in blen] 
    if flags['F']:
        nads=["[&AD={:#5.3f},ADS={:#5.3f}]:".format(x,x) for x in nad]
    elif flags['G']:
        bads=["[&AD={:#5.3f},ADS={:#5.3f}]:".format(x,x) for x in ad]
        nads=["[&AD={:#5.3f},ADS={:#5.3f}]:".format(x,x) for x in nad]
        blens1=[prec.format(x) for x in rlen] 
        blens2=[prec.format(x) for x in rlen2] 
    else:
        nads=[":"]*nnode
    clen=0
    while len(nstck)>1:
        k=nstck[-1]
        #gdbug(3,k,nwk,nwk[k])
        if nstck.count(k)>1:
            #gdbug(9,nstck)
            madlog(bug)
        if k<notu:
            lens=blens[k]            
            if flags['G']:
                if rlen[k]==0:
                    nwk[k]="{}{}{}".format(labs[k],nads[k],blens2[k])
                elif rlen2[k]==0:
                    nwk[k]="{}{}{}".format(labs[k],nads[k],blens1[k])
                else:
                    nwk[k]="({}{}{}){}{}".format(labs[k],nads[k],blens1[k],bads[k],blens2[k])
            else:            
                nwk[k]="{}{}{}".format(labs[k],nads[k],blens[k])
            nstck.pop()
            #gdbug(3,'MM',k,nwk,blens,blen,nstck)
            clen+=blen[k]
            continue
        elif len(nwk[k])==1:  
            for b in [0,1,2]:
                if nodes[k,b]==nwk[k][0]:
                    continue
                c=nodes[k,b]
                nwk[c]=[k]
                nstck.append(c)
                nwk[k].append(c)
            continue
        else:
            b=nwk[k]
            f,c,d=b
            nwk[k]="{},{}".format(nwk[c],nwk[d])
            nwk[c],nwk[d]=["",""]
            if trip[k,f]==2:
                k1=k
            else:
                k1=f     
            #gdbug(3,k,b[0],k1,blen[k1])
            if flags['G']:
                if blen[k1]==0:
                    pass
                else:        
                    if rlen[k1]==0:
                        nwk[k]="({}){}{}".format(nwk[k],nads[k],blens2[k1])
                    elif rlen2[k1]==0:
                        nwk[k]="({}){}{}".format(nwk[k],nads[k],blens1[k1])
                    else:
                        if k1==k:
                            nwk[k]="(({}){}{}){}{}".format(nwk[k],nads[k],blens1[k1],bads[k1],blens2[k1])
                        else:
                            nwk[k]="(({}){}{}){}{}".format(nwk[k],nads[k],blens2[k1],bads[k1],blens1[k1])
            elif blen[k1]>0 or not flags['P']:
                nwk[k]="({}){}{}".format(nwk[k],nads[k],blens[k1])
            #gdbug(3,"qq {}".format(blen[k1]))                 
            nstck.pop()
            clen+=blen[k1]
        #gdbug(1,'ww',k,nwk,blens,blen,nstck)
        #<--while stck    
    tol=(tlen-clen)/tlen
    #gdbug(1,"tol=",tol)
    newnwk="({},{})".format(nwk[nodes[rootnode,0]],nwk[nodes[rootnode,2]])
    if abs(tol)>0.0000001:          
                                    #gdbug(0,nodes,newnwk,nodes[rootnode,:],[rootnode,notu,newnwk.count(':')])
                                    madlog(bug+
                                    " tlen {} clen {} tol {}".format(tlen,clen,tol))
    #gdbug(2,newnwk,"\n----\n",nwk)
    return newnwk
#<--tree2nwk()
#========================================


#========================================
#====== fimd minimal values, generate AI, reroot
def mad_output():
    global nodes,blen,trip,labs,ofn,rlen,rlen2,ad,nad
    #--------- find mads        
    mad=min(ad)
    if mad<0.001: madlog("Warning: MAD=={:.5g} is too good to be true.".format(mad))
    roots = [i for i, x in enumerate(ad) if x <=(mad*madtol)]
    nroots=len(roots)
    if nroots==1:
        ai=mad/sorted(ad)[1] 
    else:
        ai=1.0
    gdbug(1,mad,roots,ai)
    #----- output:
    #---- detach rootnode
    a,b=list(nodes[rootnode,0:2])
    #gdbug(1,nodes,a,b,blen,rlen,rlen2)
    nodes[a,2]=b
    nodes[b,2]=a
    blen[b]=blen[a]
    rlen[b]=rlen2[a]
    rlen2[b]=rlen[a]
    ad[b]=ad[a]
    ad[rootnode]=mad
    nad[rootnode]=mad
    if nroots>1:
        saved=[deepcopy(v) for v in [labs,blen,nodes,trip,ad,nad,rlen,rlen2]]
    rooted=""
    #gdbug(1,nodes,blen,rlen,rlen2)
    for r in range(nroots):
        if flags['N'] and nroots>1:
            r=np.argmin([ccv[j] for j in roots])
        i=roots[r] 
        if i==rootnode:             madlog(bug+' (rootnode in roots.)')

        #reattach rootnode at inferred root    
        j=nodes[i,2]
        k=trip[j,i]
        #gdbug(2,[i,j,k])       
        #gdbug(3,nodes,blen,labs)
        nodes[i,2]=rootnode
        nodes[j,k]=rootnode
        nodes[rootnode,0]=i
        nodes[rootnode,1]=-1
        nodes[rootnode,2]=j
        trip[i,rootnode]=2
        trip[j,rootnode]=k
        trip[rootnode,i]=0
        trip[rootnode,j]=2
        if k==2:
            k1=j
            k2=rootnode
        else:
            k1=rootnode
            k2=j
        #blen[k1]=blen[i]-rlen[i]           #new branch length
        blen[k1]=rlen2[i]           #new branch length
        blen[i]=rlen[i]                    #new branch length 
        rlen[rootnode]=0 
        rlen2[rootnode]=rlen2[i] 
        rlen[k1]=0
        rlen2[k1]=rlen2[i]
        rlen2[i]=0
        if blen[k1]==0 or blen[i]==0: madlog("Warning: Root is polytomous.")
        #gdbug(1,nodes,blen,labs)
            #------- make newick string
        rootstr=tree2nwk().replace(at,"@")
        s="[MAD={:#5.3f}_AI={:#5.3f}_CCV={:#.3g}%_N={}/{}]".format(mad,ai,ccv[i]*100,r+1,nroots)
        madlog(">> "+s)
        if flags['S']:
            rootstr+=s 
        if flags['F'] or flags['G']:
            s='[&AD={:#5.3f},ADS="MAD={:#5.3f}",STT="CCV={:#.3g}%"]:{:.3g})[&ADS="AI={:#5.3f}"]'.format(mad,mad,ccv[i]*100,depth[i]*0.1,ai)
            rootstr='tree tree_{} = [&R] ({}{}'.format(r+1,rootstr,s);
        rooted+=rootstr+";"
        #gdbug(2,nwkstr,rooted)        
        if flags['N']:
            break
        rooted+="\n"
        if r<nroots-1:
            labs,blen,nodes,trip,ad,nad,rlen,rlen2=[deepcopy(v) for v in saved]
    #<--r loop
    #--------            
    rccv=", ".join(["{:#.3g}%".format(ccv[i]*100) for i in roots])
    sttstr="\nMinimal ancestor deviation, MAD = {:#5.3f}\n".format(mad)+ \
           "           Ambiguity index,  AI = {:#5.3f}\n".format(ai)+ \
           "                  Clock CV, CCV = {}".format(rccv)
    print(sttstr)
    s=""
    if flags['M']:
        if nroots>1: writeout(">> Tied root positions, {} rooted trees:".format(nroots))
        else:        writeout(">> Rooted tree:")
    if flags['F'] or flags['G']:
        rooted='#NEXUS\nbegin trees;\n'+rooted+"""end;
begin figtree;
	set appearance.branchColorAttribute="AD";
	set appearance.branchLineWidth=6.0;
	set colour.scheme.ad="AD:HSBContinuous{{false,false,0,1},0.5,1.0,1.0,0.9,1.0,0.75,false}";
	set legend.attribute="AD";
	set legend.isShown=true;
	set nodeLabels.displayAttribute="ADS";
	set nodeLabels.isShown=true;
	set branchLabels.displayAttribute="STT";
	set branchLabels.isShown=true;
	set rectilinearLayout.alignTipLabels=true;
	set rectilinearLayout.curvature=9571;
	set rectilinearLayout.rootLength=0;
"""
        if flags['G']:
            rooted+="""
	set nodeShape.size=6.0;	
	set nodeShape.isShown=true;
"""
        rooted+='end;'
    writeout(rooted)
    #gdbug(1,rooted)
    if nroots>1:
        if flags['N']: print("Tied root positions,\n{} rooted trees, but just one written to {}\n".format(nroots,ofn))
        else:          print("Tied root positions,\n{} rooted trees written to {}\n".format(nroots,ofn))
    else:        print("Rooted tree written to '{}'\n".format(ofn));
#<-- mad_output       
#--------------------------   

#--------------------------   
#----- Main loop
for inwk in range(nnwk):
    if flags['M'] or flags['N']:
        print("\nAnalyzing tree {} of {} from '{}'...".format(inwk+1,nnwk,fn))
        if flags['M']: writeout(">Rooting of input tree #{}:".format(inwk+1))
    else: print("\nAnalyzing file '{}'...".format(fn))
    nwkstr=mnwk[inwk]+";"
    try:
        #---
        notu,nnode,rootnode,nodes,blen,trip,labs,prec,tlen = nwk2tree(nwkstr)
        #---
        gdbug(1,notu,nnode,rootnode,prec,tlen)
        gdbug(1,nodes,blen,trip,labs)
        prec="{"+":#.{}".format(min(prec-2,9))+"g}"
        #--- merge pre-existing root splits    
        a,b=list(nodes[rootnode,0:2])
        blen[a]+=blen[b]
        blen[b]=0
        #gdbug(1,a,b,blen)   
        #---
        enotu,npair,dist,dij,n2n=pwdist_preproc()
        #gdbug(1,enotu,npair)   
        #---
        ad,nad,rlen,rlen2,ccv,depth=ancestor_deviations() 
        #gdbug(2,ad)   
        #---
        mad_output()    
    except madError:
        pass
    except SystemExit:
        fh.close()
        raise
    except:
        madlog(bug+" (Unidentified error.)\n")
#<-- inwk loop
s="    - Please cite DOI:10.1038/s41559-017-0193\n"
if flags['M']: writeout(">>"+s)
print("\n"+s)
fh.close()
exit()
#ffu if __name__ == '__main__':