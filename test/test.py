import subprocess


def copy_file(filename):
    try:
        subprocess.run('cp ' + filename + ' .', shell=True, check=True)
    except subprocess.CalledProcessError:
        print('\033[31m'+'Error'+'\033[0m'+': Failed in preparation')
        print('Make sure that ' + filename + ' exists.')
        exit()


def copy_scfoutput(foldername):
    subprocess.run('cp tc_energy_scf.dat ' + foldername + '/.', shell=True, check=True)
    subprocess.run('cp tc_wfc_scf.dat ' + foldername + '/.', shell=True, check=True)
    subprocess.run('cp tc_scfinfo.dat ' + foldername + '/.', shell=True, check=True)


def check_tc(output, output_ref, checks_totE):
    try:
        subprocess.run('./tc++', shell=True, check=True)
    except subprocess.CalledProcessError:
        print('Error: Failed to run ./tc++')
        print('Please copy tc++ and required input files into the parent directory.')
        exit()

    with open(output) as f:
        lines = f.readlines()

    with open(output_ref) as f_ref:
        lines_ref = f_ref.readlines()

    # Total energy check
    if checks_totE:
        lines_ref_totE = [line for line in lines_ref if 'Total energy =' in line]
        lines_ref_totE_split = lines_ref_totE[len(lines_ref_totE)-1].split()
        totE_ref = [lines_ref_totE_split[4]]
        lines_ref_totE_split = lines_ref_totE[len(lines_ref_totE)-2].split()
        totE_ref.append(lines_ref_totE_split[4])
        
        print('  Reference total energies (two SCF loops)', totE_ref)
        
        lines_totE = [line for line in lines if 'Total energy =' in line]
        lines_totE_split = lines_totE[len(lines_totE)-1].split()
        totE = [lines_totE_split[4]]
        lines_totE_split = lines_totE[len(lines_totE)-2].split()
        totE.append(lines_totE_split[4])

        print('  Calculated total energies (two SCF loops)', totE)
        
        if (abs(float(totE[0]) - float(totE_ref[0])) < 1e-5 and abs(float(totE[1]) - float(totE_ref[1])) < 1e-5):
            print('\033[32m'+'  Total-energy check passed.'+'\033[0m')
        else:
            print('\033[31m'+'Error'+'\033[0m'+': Total-energy difference is larger than 1e-5')
            exit()

    # Orbital energy check
    lines_ref_orbE = [line for line in lines_ref if '   6 ' in line]
    lines_ref_orbE_split = lines_ref_orbE[len(lines_ref_orbE)-1].split() # final line
    lines_ref_orbE_split2 = lines_ref_orbE[len(lines_ref_orbE)-2].split()
    orbE_ref = [lines_ref_orbE_split[1], lines_ref_orbE_split2[1]]
    print('  Reference orbital energies', orbE_ref)

    lines_orbE = [line for line in lines if '   6 ' in line]
    lines_orbE_split = lines_orbE[len(lines_orbE)-1].split() # final line
    lines_orbE_split2 = lines_orbE[len(lines_orbE)-2].split()
    orbE = [lines_orbE_split[1], lines_orbE_split2[1]]
    print('  Calculated orbital energies', orbE)

    if (abs(float(orbE[0]) - float(orbE_ref[0])) < 1e-5 and abs(float(orbE[1]) - float(orbE_ref[1])) < 1e-5):
        print('\033[32m'+'  Orbital-energy check passed.'+'\033[0m')
    else:
        print('\033[31m'+'Error'+'\033[0m'+': Orbital-energy difference is larger than 1e-5')
        exit()


print('\033[35m'+'Start test calculations (will take few minutes).'+'\033[0m')
print('\033[31m'+'Note! DO NOT use the input files (including pseudopot.) provided here for your research.')
print('These were made only for test calculation...'+'\033[0m'+'\n')

# common in the following examples
output = './output.out'

## Example 1
print('\033[35m'+'Test 1/5'+'\033[0m'+': HF SCF calculation for bulk Si')
copy_file('./example1/input.in')
check_tc(output, './example1/output.out', True)
print('\033[35m'+'Test 1/5'+'\033[32m'+' passed.'+'\033[0m'+'\n')

# copy output files for the subsequent tests...
copy_scfoutput('./example4')
copy_scfoutput('./example6')

## Example 2
print('\033[35m'+'Test 2/5'+'\033[0m'+': TC SCF calculation for bulk Si')
copy_file('./example2/input.in')
check_tc(output, './example2/output.out', True)
print('\033[35m'+'Test 2/5'+'\033[32m'+' passed.'+'\033[0m'+'\n')

## Example 3
print('\033[35m'+'Test 3/5'+'\033[0m'+': BITC SCF calculation for bulk Si')
copy_file('./example3/input.in')
check_tc(output, './example3/output.out', True)
print('\033[35m'+'Test 3/5'+'\033[32m'+' passed.'+'\033[0m'+'\n')

## Example 4
print('\033[35m'+'Test 4/5'+'\033[0m'+': HF SCF calculation for bulk Si (restarted calc.)')
copy_file('./example4/input.in')
copy_file('./example4/tc_energy_scf.dat')
copy_file('./example4/tc_wfc_scf.dat')
copy_file('./example4/tc_scfinfo.dat')
check_tc(output, './example4/output.out', True)
print('\033[35m'+'Test 4/5'+'\033[32m'+' passed.'+'\033[0m'+'\n')

## Example 6
print('\033[35m'+'Test 5/5'+'\033[0m'+': HF BAND calculation for bulk Si')
copy_file('./example6/input.in')
copy_file('./example6/tc_energy_scf.dat')
copy_file('./example6/tc_wfc_scf.dat')
copy_file('./example6/tc_scfinfo.dat')
check_tc(output, './example6/output.out', False)
print('\033[35m'+'Test 5/5'+'\033[32m'+' passed.'+'\033[0m'+'\n')

print('Remove temporary files...')
subprocess.run('rm ./tc_energy_scf.dat', shell=True)
subprocess.run('rm example*/tc_energy_scf.dat', shell=True)
subprocess.run('rm ./tc_wfc_scf.dat', shell=True)
subprocess.run('rm example*/tc_wfc_scf.dat', shell=True)
subprocess.run('rm ./tc_energy_band.dat', shell=True)
subprocess.run('rm ./tc_wfc_band.dat', shell=True)
subprocess.run('rm ./tc_scfinfo.dat', shell=True)
subprocess.run('rm example*/tc_scfinfo.dat', shell=True)
subprocess.run('rm ./tc_bandplot.dat', shell=True)
subprocess.run('rm ./input.in', shell=True)
subprocess.run('rm ./output.out', shell=True)

print('\033[35m'+'Since these tests are done in serial calculation, it is recommended to check MPI parallelization does not change the results by yourself...'+'\033[0m')

