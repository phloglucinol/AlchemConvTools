
import json
from pymol import cmd

with open('groups.json', 'r') as f:
    groups = json.load(f)['init_test_group']                                                                             
cmd.load('tleap-mol.pdb', 'ligand')
# cmd.select('test', 'id 1| id 2')
cmd.set ( 'ray_trace_mode', 3 )
cmd.set ( 'bg_rgb', 'white' )
cmd.set ( 'ray_opaque_background', 0 )
cmd.set ( 'valence', 0 )
cmd.set ( 'cartoon_fancy_helices', 1 )
cmd.set ( 'ray_shadow', 0 )
cmd.set ( 'cartoon_transparency', 0.5 )
color_dict = {'group1': 'salmon', 'group2':'deepblue', 'group3': 'forest', 'group4': 'sand', 'group5': 'pink', 'group6' :'cyan', 'group7': 'orange', 'group8': 'white', 'group9': 'yellow', 'group10': 'magenta', 'group11': 'limon', 'group12': 'red', 'group13': 'bluewhite', 'group14': 'greencyan', 'group15': 'lightblue', 'last_group': 'gray70'}
count=0
for key, value in groups.items():
   cmd.select(key, '|'.join(['id {}'.format(atom_id) for atom_id in value]))
   #cmd.color(color_dict[key], key)
   count +=1
   if count == len(groups):
     cmd.color('gray70', key)
   else:
     cmd.color(color_dict[key], key)
cmd.zoom('ligand', complete=1)
cmd.png('groups.png' , width=1693, height=1693, dpi=600,  ray=1 )





