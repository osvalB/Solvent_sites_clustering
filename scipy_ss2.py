#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 13:16:19 2019

@author: osvaldo

El siguiente script se usa para calcular sitios de solvente 
con el algoritmo descripto en https://watclust.wordpress.com/methodology/

Requiere las librerias:

	scipy, networkxy pytraj

Ejecutar como python scipy_ss2.py para comprender los argumentos requeridos

La nomeclatura del "binding site" debe estar en formato AMBER.

Atención: Para un correcto cálculo de la energía del sitio debe estar 
la densidad del cosolvente en # moleculas / (A^3)

"""
from __future__ import division

import numpy as np
import pytraj as pt
import os,time
from scipy import spatial
from itertools import product, groupby
from collections import Counter
import argparse
import networkx
from networkx.algorithms.components.connected import connected_components


#binding_site_resids    =    ":143,41,44,45,48,49,53,56,59,73,74,77,78,79,81,84,85,86,87,89,90,91,93,94"
ap = argparse.ArgumentParser()

ap.add_argument("-nc",     "--netcdf",              required=True)
ap.add_argument("-p",      "--prmtop",              required=True)
ap.add_argument("-bs",     "--binding_site_resids", required=True)
ap.add_argument("-dist",   "--dist_threshold",      required=True)
ap.add_argument("-outdir", "--outputdir",           required=True)

ap.add_argument("-r",       "--reference",      default="ref.pdb")
ap.add_argument("-solv",    "--solvent",        default="ETA")
ap.add_argument("-s",    "--step",        default=1)
ap.add_argument("-probe",   "--probe_atom",     default="C1")
ap.add_argument("-ff",      "--first_frame",    default="0")
ap.add_argument("-lf",      "--last_frame",     default="999")
ap.add_argument("-dr",      "--dr",             default="0.1")
ap.add_argument("-watNmin", "--watnumber_min",  default="100")
ap.add_argument("-pop",     "--pop",            default="0.9")
ap.add_argument("-WFRr",    "--WFRr",           default="0.6")
ap.add_argument("-dens",    "--density",        default="0.00224")

args = vars(ap.parse_args())
 
# display a friendly message to the user

print(args)

#######################################

#User defined:

binding_site_resids	=	args["binding_site_resids"]#":1-85"
ref					=	args["reference"]#"ref.pdb"
netcdf				=	args["netcdf"]#"center_md_1.nc"
parm				=	args["prmtop"]#"strip.mdm2._1.prmtop"
solvent				=	args["solvent"]#"ETA"
probe				=	args["probe_atom"]#"C1"
dist_threshold		=	float(args["dist_threshold"])#0.3
watnumber_min		=	int(args["watnumber_min"])#29
first_frame			=	int(args["first_frame"])#0
last_frame			=	int(args["last_frame"])#299
step				=	int(args["step"])

dr 					=	float(args["dr"])#0.1 
WFRr 				=	float(args["WFRr"])#0.6 
pop 				=	float(args["pop"])#0.90 
densidad 			=	float(args["density"])#0.0334 

def pdb_to_dict(pdb):
    serial = 1
    atom_res_dict    =    {}
    with open(pdb,"r") as p:
        for l in p:
            if "ATOM" in l or "HETATM" in l:
                atom_res_dict[serial]    =    l[23:26]
                serial                     +=     1
    return(atom_res_dict)

def pdb_writer(file,atomnumber,atomtype,resname,resid,x,y,z,PFP,R90):
    f = open(file,"a+")
    f.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f} {:5.1f}{:6.2f}          {:>2s}".format(
        "ATOM",atomnumber,atomtype," ",resname," ",resid," ",\
        x,y,z,PFP,R90,""))
    f.write('\n')
    f.write('TER')
    f.write('\n')
    f.close()
    
def are_near(atom1,atom2,dist):
    dist2 = dist*dist
    x1     = atom1
    x2     = atom2
    diff   = x1 - x2
    dd     = np.sum(diff**2, axis=0)
    return(dd <= dist2)    
    
ref                 =    pt.load(ref)
traj                =     pt.load(netcdf,parm)

traj                =    traj.autoimage()

frames              = [x for x in range(first_frame,last_frame+1,step)]

# Align trajectory to reference
mask                 =    "("+binding_site_resids+")&(!(@H))"
pt.align(traj,mask=mask,ref=ref,ref_mask=mask)
pt.rmsd(traj, mask=mask,ref=ref,ref_mask=mask,update_coordinate=True,nofit=False)

# Make serial ID and residue number dictionary
scratch_pdb_name    = "scratchPDB-XAAZSAASA.pdb"
pt.write_traj(scratch_pdb_name,traj,frame_indices = [0],overwrite=True)
dicc                 = pdb_to_dict(scratch_pdb_name)
os.system("rm -f " + scratch_pdb_name)

def to_graph(l):
    G = networkx.Graph()
    for part in l:
        # each sublist is a bunch of nodes
        G.add_nodes_from(part)
        # it also imlies a number of edges:
        G.add_edges_from(to_edges(part))
    return G

def to_edges(l):
    """ 
        treat `l` as a Graph and returns it's edges 
        to_edges(['a','b','c','d']) -> [(a,b), (b,c),(c,d)]
    """
    it = iter(l)
    last = next(it)

    for current in it:
        yield last, current
        last = current    

def clusterize(probe,solvent,outputdir):

    positions             =    []
    metadata             =    []
    mask_selection     =     "(@"+probe+")&(:"+solvent+")"+"&("+binding_site_resids+"<:3.8)"
    n = 0
    
    start = time.time()

    for frame in frames:
        n                 +=     1
        traj.top.set_reference(traj[frame])
        top             =     traj.top.select(mask_selection)
        xyz                =    traj[frame,mask_selection].xyz
    
        for a in range(len(top)):
            x,y,z                 = float(xyz[a][0]),float(xyz[a][1]),float(xyz[a][2])
            serial_i    =    top[a]+1
            resid             =    int(dicc[serial_i])
            positions.append((x,y,z))
            metadata.append([probe,frame,solvent,resid,serial_i])    
   
    os.system("mkdir -p " + outputdir)
    os.chdir(outputdir)

    tempf = open("params.dat","w")
    for (i,k) in args.items():
        tempf.write("{},{}\n".format(i,k))
    tempf.close()
 
    points = np.asarray(positions)
    
    data_time = time.time() - start

    tree = spatial.cKDTree(points,leafsize=50)
    
    start2 = time.time()

    l = [(tree.query_ball_point(p,dist_threshold)) for p in points]
    
    ##### Primer poda ##########

    lists = [x for x in l if len(x) > 1]

    query_time = time.time() - start2

    join_time = time.time()

    G = to_graph(lists)
    l = list(connected_components(G))   
 
    join_time = time.time() - join_time

    clusters = [list(x) for x in l if len(x) >= watnumber_min]
    
    file = open("watclust.dat","w")
    file.write("WS\t#TotWat\tWFP(r={})\tr(pop={})\t#Aguas_distintas\n".format(WFRr,pop))  
               
    cluster_n = 0   
            
    for cl in clusters:
        posT = points[cl]
        ags = []
        for cx in cl:
            if metadata[cx][3] not in ags:
                ags.append(metadata[cx][3])

        ag_dist = len(ags)

        centroid = posT.mean(axis=0)
        selWFR = [at for at in posT if are_near(at,centroid,WFRr)]
        WFR = len(selWFR) / (n * 4/3 * 3.141593 * WFRr * WFRr * WFRr * densidad)
        cluster_n += 1
        name = "cluster_"  + str(cluster_n) +".pdb"
        os.system("rm -f " + str(name))
    
        for pos,data in zip(posT,[metadata[i] for i in cl]):
            pdb_writer(name,data[4],probe,data[2],\
                       data[3],pos[0],pos[1],pos[2],0,0)		
    
        popr		=	0
        r 			=	0
        os.system("rm -f " + "watr"+str(cluster_n)+".dat")
        file = open("watr"+str(cluster_n)+".dat","w")
    
        while popr < 1:
            r 		=	r + dr
            if r == dr:
                selgr 	=	[atom for atom in posT if are_near(atom,centroid,r)]
            else:
                selgr 	=	[atom for atom in posT if (are_near(atom,centroid,r) and not are_near(atom,centroid,(r-dr)))]
                selwfrr 	=	[atom for atom in posT if are_near(atom,centroid,r)]
    	  			
                gr 			=	len(selgr) / (4/3 * 3.141593 * ( (r * r * r) - ((r-dr) * (r-dr) * (r-dr)) ) * densidad * n)
    				
                wfr 		=	len(selwfrr) / (n * 4/3 * 3.141593 * (r * r * r) * densidad)
    
                popr 		=	len(selwfrr) * 1 / len(cl)
    
                file.write("{} {} {}\n".format((r - dr) / 2,gr,wfr))
    
            if popr <= pop:
                rpop 	=	r
    
        file.close()
    
        pdb_writer("watcent.pdb",cluster_n,probe,metadata[0][2],\
                   cluster_n,centroid[0],centroid[1],centroid[2],0,0)
    
        file = open("watclust.dat","a+")
        file.write("WS{}\t{}\t{}\t{}\t{}\n".format(cluster_n,len(cl),WFR,round(rpop,2),ag_dist))
        file.close()
    print('It took '+str(time.time()-start)+' seconds.')

    print("{},{},{},{},{}".format(len(points),time.time()-start,data_time,query_time,join_time))

print("")

clusterize(probe,solvent,args["outputdir"])#"C1")    
    
    
    
