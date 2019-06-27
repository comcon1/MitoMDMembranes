# Copyright 2019 Simon & Garfunkel
import colorsys,sys
from pymol import cmd
import numpy as np

def rithm():
    cmd.bg_color("white")
#    cmd.load("view*")
    cmd.set('defer_builds_mode', 3)
    cmd.set("antialias", 5)
    cmd.set("ray_trace_mode", 6)
    cmd.hide()
    cmd.select("mito_lipids", "resname POPC or resname SLPE or resname PLPC")
    cmd.select("tlx", "resname TLX*")
    cmd.color("hotpink", "resname 3PP")
    cmd.color("deeppurple", "resname BDP")
    cmd.color("marine", "resname BTL")
    cmd.show_as("sticks", "mito_lipids or tlx")
    cmd.select("lipids_H", "hydro and (resname POPC or resname SLPE or resname PLPC or resname TLX*)")
    cmd.hide("everything", "lipids_H")
    cmd.select("bodipy", "resname BDP or resname 3PP or resname BTL")
    cmd.show_as("sticks", "bodipy")
    cmd.set("stick_radius", 0.6, "bodipy")
#  cmd.show("spheres", "bodipy")
#  cmd.set("sphere_transparency",0.2)
#    cmd.set("stick_color", "white", "bodipy")
    cmd.select("perox_o", "name OP* and tlx")
    cmd.select("oxy", "name O* and (mito_lipids or tlx)")
    cmd.select("phos", "name P* and (mito_lipids or tlx)")
    cmd.select("carbon_MITO", "name *C* and mito_lipids")
    cmd.select("carbon_TLX", "name *C* and tlx")
    cmd.color("gray80", "carbon_MITO")   
    cmd.color("tv_green", "carbon_TLX")  
    cmd.select("water", "resname SOL") 
    cmd.select("na", "resname NA") 
    cmd.select("cl", "resname CL") 
    cmd.color("salmon", "cl") 
    cmd.color("yellow", "na") 
    cmd.show_as("spheres", "na or cl")
    cmd.set("sphere_scale", 0.5, "na or cl")
    cmd.show_as("surface", "water")
    cmd.set("transparency", 0.5)
    cmd.color("aquamarine", "water")
    cmd.set("ray_shadows", "off")
    cmd.set("surface_solvent", "on")
    cmd.set("transparency_mode", 1)

#    run drawBoundingBox.py
#drawBoundingBox padding=-13.0, linewidth=2.0, r=0, g=0, b=0    
