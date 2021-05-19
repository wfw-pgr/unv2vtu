import os, sys
import numpy                          as np
import nkVTKRoutines.construct__uGrid as cug

# ========================================================= #
# ===  unv2vtu.py  ( ideas-universal format to vtu )    === #
# ========================================================= #

def unv2vtu( mshFile=None, fldFiles=None, vtkFile=None, elmFile=None ):

    # ------------------------------------------------- #
    # --- [1] input files                           --- #
    # ------------------------------------------------- #
    
    if ( mshFile  is None ): sys.exit( "[unv2vtu.py] mshFile  == ??? [ERROR] " )
    if ( fldFiles is None ): sys.exit( "[unv2vtu.py] flsFiles == ??? [ERROR] " )
    if ( elmFile  is None ): sys.exit( "[unv2vtu.py] elmFile  == ??? [ERROR] " )
    if ( vtkFile  is None ): vtkFile = "dat/out.vtu"
    
    # ------------------------------------------------- #
    # --- [2] data loading                          --- #
    # ------------------------------------------------- #
    #  -- [2-1]  mesh loading                       --  #
    mesh     = load__meshFile( mshFile=mshFile )
    physNums = load__elemFile( inpFile=elmFile )
    elems    = mesh["mesh"] - 1
    nodes    = mesh["node"]
    index    = mesh["mesh_info"][:,0] - 1
    sortidx  = np.argsort( index )
    index    = index[sortidx]
    elems    = elems[sortidx]
    physNums = physNums[sortidx]
    nElems   = elems.shape[0]
    # -- node number begin from 1 not 0 -- #
    
    #  -- [2-2]  fields loading                     --  #
    fields   = []
    for fldFile in fldFiles:
        fields.append( load__fieldFile( fldFile=fldFile ) )

    #  -- [2-3]  decomposit field into numpy array  --  #
    fkeys    = []
    fieldL   = []
    indexL   = []
    for field in fields:
        hkeys  = list( field["field"].keys() )
        fkeys  = fkeys  + hkeys
        for ik,key in enumerate( hkeys ):
            fieldL.append( field["field"][key] )
            indexL.append( np.searchsorted( index, field["index"]-1 ) )
            print( "[unv2vtu.py] key :: {0:10} is added... ".format( key ) )
            # -- index order should begin 0 -- #
    fkeys .append( "physNum" )
    fieldL.append( physNums[:,1] )
    indexL.append( physNums[:,0] )
    nField   = len( fkeys )

    #  -- [2-4] summarize the field                  --  #
    field    = np.zeros( (nElems,nField,), dtype=np.float )
    for ik,hfield in enumerate( fieldL ):
        data             = np.zeros( (nElems,) )
        data[indexL[ik]] = hfield
        field[:,ik]      = np.copy( data )
    
    print()
    print( "[unv2vtu.py] shape of elements. == {0}".format( elems.shape ) )
    print( "[unv2vtu.py] shape of nodes.    == {0}".format( nodes.shape ) )
    print( "[unv2vtu.py] shape of field     == {0}".format( field.shape ) )
    print( "[unv2vtu.py] keys  of field     == {0}".format( fkeys       ) )
    print()
    
    # ------------------------------------------------- #
    # --- [3] reconstruct as vtu ugrid              --- #
    # ------------------------------------------------- #
    ret = cug.construct__uGrid( elems=elems, nodes=nodes, cellData=field, cellDataName=fkeys, vtkFile=vtkFile )

    # ------------------------------------------------- #
    # --- [4] write separately depending on physNum --- #
    # ------------------------------------------------- #
    physNums_list = set( physNums[:,1] )
    for phys in physNums_list:
        index = np.where( physNums[:,1] == phys )
        helem = elems[index]
        hcell = field[index]
        hnode = np.copy( nodes )
        fname = vtkFile.replace( ".vtu", "_{0}.vtu".format( phys ) )
        ret   = cug.construct__uGrid( elems=helem, nodes=hnode, cellData=hcell, cellDataName=fkeys, \
                                      vtkFile=fname )
    return()


# ========================================================= #
# ===  load__fieldFile                                  === #
# ========================================================= #

def load__fieldFile( fldFile=None ):
    
    if ( fldFile is None ): sys.exit( "[load__fieldFile] fldFile == ??" )

    # ------------------------------------------------- #
    # --- [1] load field's header                   --- #
    # ------------------------------------------------- #

    with open( fldFile, "r" ) as f:
        lines = f.readlines()
    nLines = len( lines )
    stack  = []
    for il,line in enumerate( lines ):
        nWord = len( ( line.strip() ).split() )
        word  = line.strip()
        if ( ( nWord == 1 ) and ( word == "-1" ) ):
            stack.append( il )

    if ( len( stack ) <= 0 ): sys.exit( "[load__fieldFile] can't find [ -1 / -1] lines.... ERROR " )
    stack  = np.reshape( np.array( stack ), (-1,2) )
    nBlock = stack.shape[0]

    for ik in range( nBlock ):
        ito    = stack[ik,1] - 1
        ifrom  = stack[ik,0] + 2
        fieldl = lines[ifrom:ito+1]
        items  = ( ( ( ( fieldl[0] ).split( "-" ) )[1] ).split( "," ) )
        items  = [ item.strip() for item in items ]
        index  = np.array( [ ( line.split() )[0] for line in fieldl[8::2] ], dtype=np.int   )
        field  = np.array( [ ( line.split() )    for line in fieldl[9::2] ], dtype=np.float )

    fDict = { "{0}".format( items[ik] ):field[:,ik] for ik in range( len(items) ) }
        
    return( { "index":index,"field":fDict } )


        
# ========================================================= #
# ===  load__meshFile                                   === #
# ========================================================= #

def load__meshFile( mshFile=None ):

    if ( mshFile is None ): sys.exit( "[load__meshFile] mshFile == ???" )

    # ------------------------------------------------- #
    # --- [1] load structure                        --- #
    # ------------------------------------------------- #

    with open( mshFile, "r" ) as f:
        lines = f.readlines()
    nLines = len( lines )
    stack  = []
    for il,line in enumerate( lines ):
        nWord = len( ( line.strip() ).split() )
        word  = line.strip()
        if ( ( nWord == 1 ) and ( word == "-1" ) ):
            stack.append( il )

    if ( len( stack ) <= 0 ): sys.exit( "[load__meshFile] can't find [ -1 / -1] lines.... ERROR " )
    stack  = np.reshape( np.array( stack ), (-1,2) )
    nBlock = stack.shape[0]

    
    # ------------------------------------------------- #
    # --- [2] block type                            --- #
    # ------------------------------------------------- #
    btype  = []
    for ik in range( nBlock ):
        ctype = lines[ stack[ik,0] + 1 ]
        if   ( ctype.strip() == "2411" ):
            btype.append( "node" )
        elif ( ctype.strip() == "2412" ):
            btype.append( "mesh" )

    stack  = np.copy( stack[:2,:] )
    nBlock = stack.shape[0]
    btype  = btype[:2]

    
    # ------------------------------------------------- #
    # --- [3] data fetch                            --- #
    # ------------------------------------------------- #
    for ik,bt in enumerate( btype ):

        skip_header = stack[ik,0] + 2
        skip_footer = nLines - stack[ik,1]
        
        if   ( bt == "node" ):
            ito       = stack[ik,1] - 1
            ifrom     = stack[ik,0] + 2
            nNode     = ( ito - ifrom + 1 ) // 2
            node_info = []
            node      = []
            for iN in range( nNode ):
                line1  = ( ( lines[2*iN+ifrom]   ).strip() ).split()
                line2  = ( ( lines[2*iN+ifrom+1] ).strip() ).split()
                node_info.append( [   int( sl ) for sl in line1 ] )
                node.append     ( [ float( sl ) for sl in line2 ] )
            node_info = np.array( node_info )
            node      = np.array( node     )
            
        elif ( bt == "mesh" ):

            ito       = stack[ik,1] - 1
            ifrom     = stack[ik,0] + 2
            nMesh     = ( ito - ifrom + 1 ) // 2
            mesh_info = []
            mesh      = []
            for iM in range( nMesh ):
                line1  = ( ( lines[2*iM+ifrom]   ).strip() ).split()
                line2  = ( ( lines[2*iM+ifrom+1] ).strip() ).split()
                mesh_info.append( [   int( sl ) for sl in line1 ] )
                mesh.append     ( [   int( sl ) for sl in line2 ] )
            mesh_info = np.array( mesh_info )
            mesh      = np.array( mesh      )
            
    print( "[load__meshFile] nNode == {0}, nMesh == {1}".format( nNode, nMesh ) )
    return( { "node":node, "mesh":mesh, "mesh_info":mesh_info } )



# ========================================================= #
# ===  load__elemFile                                   === #
# ========================================================= #

def load__elemFile( inpFile="dat/elem" ):

    # ------------------------------------------------- #
    # --- [1] load elem file                        --- #
    # ------------------------------------------------- #
    with open( inpFile, "r" ) as f:
        Data = np.loadtxt( f, skiprows=1 )
    Data[:,0] = Data[:,0] - 1 # 0 : elemNums
    ret = np.array( Data[:,0:2], dtype=np.int64 )
    return( ret )



# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):

    mshFile  = "dat/post_geom.unv"
    elmFile  = "dat/elem"
    fldFiles = [ "dat/magnetic.unv", "dat/magnetization.unv" ]
    vtkFile  = "dat/out.vtu"
    
    unv2vtu( mshFile=mshFile, fldFiles=fldFiles, vtkFile=vtkFile, elmFile=elmFile )
