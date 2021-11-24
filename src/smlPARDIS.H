      integer::minproc,maxproc,maxcomm,maxnodelem,maxdim,maxsv
      integer::maxnp,maxeln,maxstatev,msv,maxgrp,maxprops,maxel
      integer::maxnode,maxnset,mcn,maxns,maxdf,ndofelx
      integer::mdofe,mdofx,mcstr,mbound1,mbound2,mbound3
      integer::maxcrd,mnelx,lr,lw,nnodex,nslip,mfac,mdofproc
      integer::maxcut,maxelset,maxels,maxpress,mphase
      integer::ibrickortet,maxnrhs
      integer::mpiinfsizemax, errorcalornot
      ! ------- added by deniz
      real(8)::meshRefinErrorThrshld
      integer::maxneighbors, maxPatchElems, maxPatchNodes
      integer::MAX_HISTORY_N

      parameter(meshRefinErrorThrshld=0.80D0)
      parameter(maxneighbors=24)
      parameter(maxPatchElems=200,maxPatchNodes=300)
      parameter(MAX_HISTORY_N=4001)
      ! ----------------------------------------------
        
      parameter(minproc=1,maxproc=24)
      parameter(maxcomm=maxproc*maxproc,maxnodelem=120)
      parameter(maxdim=3,ndofelx=maxdim)
      parameter(maxnp=1,maxeln=4)
      parameter(maxstatev=386,msv=maxstatev*maxnp)
      parameter(maxgrp=1,maxprops=376,maxel=34993)
      parameter(maxsv=msv*maxel/minproc)
      parameter(maxnode=6860,maxnset=20,mcn=1)
      parameter(maxns=400)
      parameter(maxdf=ndofelx,mdofe=maxdf*maxeln)
      parameter(mdofx=ndofelx*maxnode,mcstr=mdofx)
      parameter(mbound1=10000,mbound2=10,mbound3=10)
      parameter(maxcrd=maxdim*maxnode,mnelx=maxeln*maxel)
      parameter(lr=8,lw=9)
      parameter(nnodex=maxnode,nslip=100)
      parameter(mfac=8,mdofproc=mdofx*mfac/minproc)
        parameter(maxcut=20)
      parameter(maxelset=24,maxels=2000,maxpress=24,mphase=2)
      parameter(ibrickortet=2)
        parameter(mpiinfsizemax=4200)

