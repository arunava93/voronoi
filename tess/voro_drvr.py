from .voro_helper import Timer

def wrt_main(voro,no_state):
    
    lines=[]
    lines.append('INCLUDE ../models/getstressNph.f90')
    lines.append('STATE_VAR '+str(voro.n*(no_state+6)+2))
    lines.append('LOAD_TYPE 4')
    
    return lines

def wrt_materials(voro,no_state,mat,lines):
    
    lines.append('BEGIN_MATERIAL')
    lines.append(str(float(voro.n))+' '+str(float(len(mat)+3))+' '+
                 str(float(no_state))+' 6.0')
    for i in range(len(voro.gids)):        
        mat_tot=voro.euler[i]+mat[:]
        lines.append(' '.join('{:.5e}'.format(m) for m in mat_tot))
    lines.append(' '.join('{:.5e}'.format(v/voro.v) for v in voro.vols))
    lines.append('END_MATERIAL')

    return lines

def wrt_steps(voro,D,inc_size,lines):
    
    lines.append('BEGIN_LOAD')

    ss=[str(i+1) for i in range(len(D)) if D[i]]
    lines.append(' '.join(str(s) for s in ss))

    sv=' '.join(str(a)[1:-1] if a else '0.0' for a in D)
    lines.append(str(int(1/inc_size))+' '+sv+' 1.0d0')
    
    lines.append('END_LOAD')
    
    return lines

def driver_export(voro,driver_input):
    
    no_state=driver_input.get('no_state_vars',0)
    mat=driver_input.get('mat',[])
    D=driver_input.get('D',[])
    inc_size=driver_input.get('inc_size',0.01)

    t=Timer()
    t.start()

    lines=wrt_main(voro,no_state)
    lines=wrt_materials(voro,no_state,mat,lines)
    lines=wrt_steps(voro,D,inc_size,lines)
    
    with open(voro.fname+".inp","w") as f:
        f.write('\n'.join(lines))
        f.write('\n')

    t.stop('Elapsed time in Abaqus export')