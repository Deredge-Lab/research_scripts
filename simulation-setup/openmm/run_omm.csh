#!/bin/csh

set init = step3_input
set min1_prefix = min_wat
set min2_prefix = min_prot
set equi1_prefix = heat_NVT
set equi2_prefix = eq1_NPT
set equi3_prefix = eq2_NPT
set equi4_prefix = eq3_NPT
set equi5_prefix = eq4_NPT
set equi6_prefix = eq5_NPT
set prod_prefix = prod
set prod_step   = prod

### Minimization stage 1 (relax solvent)
set input_param = "-t toppar.str -p ${init}.psf -c ${init}.crd  -b sysinfo.dat"
python -u openmm_run.py -i ${min1_prefix}.inp ${input_param} -orst ${min1_prefix}.rst > ${min1_prefix}.out

### Minimization stage 2 (relax protein)
set input_param = "-t toppar.str -p ${init}.psf -c ${init}.crd -irst ${min1_prefix}.rst"
python -u openmm_run.py -i ${min2_prefix}.inp ${input_param} -orst ${min2_prefix}.rst > ${min2_prefix}.out

### NVT heating/equilibration stage 1 (protein heavy atom restraints = 10 kcal/mol/A**2) 
set input_param = "-t toppar.str -p ${init}.psf -c ${init}.crd -irst ${min2_prefix}.rst"
python -u openmm_run.py -i ${equi1_prefix}.inp ${input_param} -orst ${equi1_prefix}.rst -odcd ${equi1_prefix}.dcd > ${equi1_prefix}.out

### NPT equilibration stage 2 (restraints = 10.0) 
set input_param = "-t toppar.str -p ${init}.psf -c ${init}.crd -irst ${equi1_prefix}.rst"
python -u openmm_run.py -i ${equi2_prefix}.inp ${input_param} -orst ${equi2_prefix}.rst -odcd ${equi2_prefix}.dcd > ${equi2_prefix}.out

### NPT equilibration stage 3 (restraints = 5.0) 
set input_param = "-t toppar.str -p ${init}.psf -c ${init}.crd -irst ${equi2_prefix}.rst"
python -u openmm_run.py -i ${equi3_prefix}.inp ${input_param} -orst ${equi3_prefix}.rst -odcd ${equi3_prefix}.dcd > ${equi3_prefix}.out

### NPT equilibration stage 4 (restraints = 2.5) 
set input_param = "-t toppar.str -p ${init}.psf -c ${init}.crd -irst ${equi3_prefix}.rst"
python -u openmm_run.py -i ${equi4_prefix}.inp ${input_param} -orst ${equi4_prefix}.rst -odcd ${equi4_prefix}.dcd > ${equi4_prefix}.out

### NPT equilibration stage 5 (restraints = 1.0) 
set input_param = "-t toppar.str -p ${init}.psf -c ${init}.crd -irst ${equi4_prefix}.rst"
python -u openmm_run.py -i ${equi5_prefix}.inp ${input_param} -orst ${equi5_prefix}.rst -odcd ${equi5_prefix}.dcd > ${equi5_prefix}.out

### NPT equilibration stage 6 (restraints = 0.1)
set input_param = "-t toppar.str -p ${init}.psf -c ${init}.crd -irst ${equi5_prefix}.rst"
python -u openmm_run.py -i ${equi6_prefix}.inp ${input_param} -orst ${equi6_prefix}.rst -odcd ${equi6_prefix}.dcd > ${equi6_prefix}.out

### Production (no restraints)
### Total production time in ns = cntmax*[(nsteps*dt)/1000]
### See prod.inp for nsteps and dt; see below for cntmax
set cnt = 1
set cntmax = 10

while ( ${cnt} <= ${cntmax} )
    @ pcnt = ${cnt} - 1
    set istep = ${prod_step}_${cnt}
    set pstep = ${prod_step}_${pcnt}
    if ( ${cnt} == 1 ) set pstep = ${equi6_prefix}
    set input_param = "-t toppar.str -p ${init}.psf -c ${init}.crd -irst ${pstep}.rst"
    python -u openmm_run.py -i ${prod_prefix}.inp ${input_param} -orst ${istep}.rst -odcd ${istep}.dcd > ${istep}.out
   @ cnt += 1
end


