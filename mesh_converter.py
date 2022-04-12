import numpy as np
import os
import bempp.api

def convert_pqr2xyzr(mesh_pqr_path, mesh_xyzr_path):
    pqr_file = open(mesh_pqr_path, 'r')
    pqr_data = pqr_file.read().split('\n')
    xyzr_file = open(mesh_xyzr_path, 'w')
    for line in pqr_data:
        line = line.split()
        if len(line) == 0 or line[0] != 'ATOM':
            continue
        xyzr_file.write(line[5]+"\t"+line[6]+"\t"+line[7]+"\t"+line[9]+"\n")
    pqr_file.close()
    xyzr_file.close()
    
def generate_msms_mesh(mesh_xyzr_path, output_dir, output_name, density, probe_radius):
    path = os.path.join(output_dir, output_name)
    command = "msms -if "+mesh_xyzr_path+" -of "+path+" -p "+str(probe_radius)+" -d "+str(density)+" -no_header"
    os.system(command)
    
def import_msms_mesh(mesh_face_path, mesh_vert_path):
    face = open(mesh_face_path, 'r').read()
    vert = open(mesh_vert_path, 'r').read()

    faces = np.vstack(np.char.split(face.split('\n')[0:-1]))[:, :3].astype(int) - 1
    verts = np.vstack(np.char.split(vert.split('\n')[0:-1]))[:, :3].astype(float)

    #grid = bempp.api.grid_from_element_data(verts.transpose(), faces.transpose())
    grid = bempp.api.Grid(verts.transpose(), faces.transpose())
    return grid

def convert_msms2off(mesh_face_path, mesh_vert_path, mesh_off_path):
    face = open(mesh_face_path, 'r').read()
    vert = open(mesh_vert_path, 'r').read()

    faces = np.vstack(np.char.split(face.split('\n')[0:-1]))[:, :3].astype(int) - 1
    verts = np.vstack(np.char.split(vert.split('\n')[0:-1]))[:, :3].astype(float)

    data = open(mesh_off_path, 'w')
    data.write("OFF"+"\n")
    data.write(str(verts.shape[0])+" "+str(faces.shape[0])+" "+str(0)+"\n")
    for vert in verts:
        data.write(str(vert[0])+" "+str(vert[1])+" "+str(vert[2])+"\n")
    for face in faces:
        data.write("3"+" "+str(face[0])+" "+str(face[1])+" "+str(face[2])+"\n")

def import_off_mesh(mesh_off_path):
    grid = bempp.api.import_grid(mesh_off_path)
    return grid

def generate_nanoshaper_mesh(mesh_xyzr_path, output_dir, output_name_temp, output_name, density, probe_radius, save_mesh_build_files):
    
    nanoshaper_dir = "C:\\APBS-3.0.0\\bin"
    nanoshaper_temp_dir = os.path.join(output_dir, "nano\\")
    mesh_dir = output_dir

    if not os.path.exists(nanoshaper_temp_dir):
        os.makedirs(nanoshaper_temp_dir)

    # Execute NanoShaper
    config_template_file = open(nanoshaper_dir+'\\config', 'r')
    config_file = open(nanoshaper_temp_dir + 'surfaceConfiguration.prm', 'w')
    for line in config_template_file:
        if 'XYZR_FileName' in line:
            path = os.path.join(mesh_dir, output_name_temp +'.xyzr')
            line = 'XYZR_FileName = ' + path + ' \n'
        elif 'Grid_scale' in line:
            line = 'Grid_scale = {:04.1f} \n'.format(density)
        elif 'Probe_Radius' in line:
            line = 'Probe_Radius = {:03.1f} \n'.format(probe_radius)

        config_file.write(line)

    config_file.close()
    config_template_file.close()

    
    os.chdir(nanoshaper_temp_dir)
    os.system(nanoshaper_dir+"\\Nanoshaper")
    os.chdir('..')
    
    os.system('move ' + nanoshaper_temp_dir + 'triangulatedSurf.vert  ' + output_name + '.vert')
    os.system('move ' + nanoshaper_temp_dir + 'triangulatedSurf.face  ' + output_name + '.face')
    
    vert_file = open(output_name + '.vert', 'r')
    vert = vert_file.readlines()
    vert_file.close()
    face_file = open(output_name + '.face', 'r')
    face = face_file.readlines()
    face_file.close()
    
    os.remove(output_name + '.vert')
    os.remove(output_name + '.face')  

    vert_file = open(output_name + '.vert', 'w')
    vert_file.write(''.join(vert[3:]))
    vert_file.close()
    face_file = open(output_name + '.face', 'w')
    face_file.write(''.join(face[3:]))
    face_file.close()

    os.chdir(nanoshaper_temp_dir)
    os.remove('exposed.xyz')
    os.remove('exposedIndices.txt')
    os.remove('stderror.txt')
    os.remove('surfaceConfiguration.prm')
    os.chdir('..')
    os.rmdir('nano')


def pqrtomesh(directory,protein,forcefield,density,probe_radius,build_mesh='yes'):
    dir_prot = directory + '\\pqr_files\\' + protein
    pf = protein +'_' + forcefield
    if density < 10.0:
        pfd = protein +'_' + forcefield + '_' +'d'+str(density)[::2]
    else:
        pfd = protein +'_' + forcefield + '_' +'d'+str(density)[:2]+'0'
    if build_mesh=='yes':
        convert_pqr2xyzr('{}/{}.pqr'.format(dir_prot,pf),'{}/{}.xyzr'.format(dir_prot,pf))
        generate_nanoshaper_mesh('{}/{}.xyzr'.format(dir_prot,pf),dir_prot,pf,pfd,density,probe_radius,False)
        grid = import_msms_mesh('{}/{}.face'.format(dir_prot,pfd),'{}/{}.vert'.format(dir_prot,pfd))
    else:
        grid = import_msms_mesh('{}/{}.face'.format(dir_prot,pfd),'{}/{}.vert'.format(dir_prot,pfd))
        
    # Read charges and coordinates from the .pqr file
    q, x_q = np.array([]), np.empty((0,3))
    molecule_file = open('{}/{}.pqr'.format(dir_prot,pf), 'r').read().split('\n')
    for line in molecule_file:
        line = line.split()
        if len(line)==0 or line[0]!='ATOM': continue
        q = np.append( q, float(line[8]))
        x_q = np.vstack(( x_q, np.array(line[5:8]).astype(float) ))
    
    return grid, q, x_q