from tess import Voro

voro_input={
  'voro_name':'voro',
  'dim':'3d',
  'no_grains':10,
  'lloyd':100,
  'size_dst':0.0,
  # 'cuts':[0,0,1,0],
  'size':[10.0,10.0,1.0],
  'periodicity':[1,1,1]
  }

mesh_input={
  'mesh_size':3.0,
  'element_order':1,
  'element_type':'tet',
  'show_ui':True,
  'save_mesh':False
  }

abaqus_input={
  'mat':[12.0,72000.0,0.0,0.3,2.86e-7,1.0e7,1.0e9,0.4,18.0],
  'no_state_vars':50,
  'non_linear':True,
  'inc_size':0.02,
  'D':[[0.1],[],[],[],[],[]],
  'L':[[],[],[],[],[],[]]
  }

voro=Voro(voro_input)
voro.plot()

voro.create_mesh(mesh_input)
voro.abaqus_export(abaqus_input)

voro=Voro('voro_name') # To read saved tesselation from file