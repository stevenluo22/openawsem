#!/usr/bin/env python3
import os
import sys
import time
import argparse
import importlib.util

from openawsem import *
from openawsem.helperFunctions.myFunctions import *

do = os.system
cd = os.chdir

def run(args):
    simulation_platform = args.platform
    platform = Platform.getPlatformByName(simulation_platform)
    if simulation_platform == "CPU":
        if args.thread != -1:
            platform.setPropertyDefaultValue("Threads", str(args.thread))
        print(f"{simulation_platform}: {platform.getPropertyDefaultValue('Threads')} threads")
    elif simulation_platform=="OpenCL":
        platform.setPropertyDefaultValue('OpenCLPlatformIndex', '0')
        platform.setPropertyDefaultValue('DeviceIndex', str(args.device))
    elif simulation_platform=="CUDA":
        platform.setPropertyDefaultValue('DeviceIndex', str(args.device))

    # if mm_run.py is not at the same location of your setup folder.
    setupFolderPath = os.path.dirname(args.protein)
    setupFolderPath = "." if setupFolderPath == "" else setupFolderPath
    proteinName = pdb_id = os.path.basename(args.protein)

    pwd = os.getcwd()
    toPath = os.path.abspath(args.to)
    checkPointPath = None if args.fromCheckPoint is None else os.path.abspath(args.fromCheckPoint)
    forceSetupFile = None if args.forces is None else os.path.abspath(args.forces)
    parametersLocation = "." if args.parameters is None else os.path.abspath(args.parameters)
    os.chdir(setupFolderPath)



    # chain=args.chain
    chain=args.chain
    pdb = f"{pdb_id}.pdb"

    if chain == "-1":
        chain = getAllChains("crystal_structure.pdb")
        print("Chains to simulate: ", chain)


    if args.to != "./":
        # os.system(f"mkdir -p {args.to}")
        os.makedirs(toPath, exist_ok=True)
        os.system(f"cp {forceSetupFile} {toPath}/forces_setup.py")
        os.system(f"cp crystal_structure.fasta {toPath}/")
        os.system(f"cp crystal_structure.pdb {toPath}/")
        # os.system(f"cp {pdb} {args.to}/{pdb}")
        # pdb = os.path.join(args.to, pdb)

    if args.fromOpenMMPDB:
        input_pdb_filename = proteinName
        seq=read_fasta("crystal_structure.fasta")
        print(f"Using Seq:\n{seq}")
    else:
        suffix = '-openmmawsem.pdb'
        if pdb_id[-len(suffix):] == suffix:
            input_pdb_filename = pdb_id
        else:
            input_pdb_filename = f"{pdb_id}-openmmawsem.pdb"
        seq=None

    if args.fasta == "":
        seq = None
    else:
        seq = seq=read_fasta(args.fasta)
        print(f"Using Seq:\n{seq}")
    # start simulation
    collision_rate = 5.0 / picoseconds
    checkpoint_file = "checkpnt.chk"
    checkpoint_reporter_frequency = 10000


    ## Configure number of timesteps and number of frames/recording interval
    # assign number of timesteps
    if args.steps < 1:
        logging.warning("--steps must be a positive integer. Reverting to default 1e7")
        total_steps = 1e7
    else:
        total_steps = args.steps
    # merge input from deprecated option --reportFrequency to option --reportInterval and set number of frames
    if args.reportFrequency != args.reportInterval:
        if args.numFrames != -1:
            logging.warning("     --reportFrequency/--reportInterval take priority over --numFrames. Ignoring user-specified --numFrames.")
        if args.reportInterval == 1000: # we've specified report frequency but not report interval
            logging.warning("     Deprecation Warning: --reportFrequency is deprecated in favor of --reportInverval, which does the same thing")
            logging.warning("     Assigning value of --reportFrequency to --reportInterval")
            report_interval = args.reportFrequency
        elif args.reportFrequency == 1000: # we've specified report interval but not report frequency
            report_interval = args.reportInterval
        else: # neither are the default value of 1000, so we interpret this as a contradiction
            logging.warning("     Deprecation Warning: --reportFrequency is deprecated in favor of --reportInverval, which does the same thing")
            logging.warning("     Reporting interval from --reportFrequency and --reportInterval disagree. Using value from --reportInterval")
            report_interval = args.reportInterval
        num_frames = int(total_steps/report_interval)
    elif args.reportFrequency == args.reportInterval == 1000: # both equal to default, so try to get number of frames
        if args.numFrames == -1: # number of frames not specified, so we can revert to default interval of 1000
            report_interval = 1000
            num_frames = int(total_steps/report_interval)
        else: # use given number of frames
            num_frames = args.numFrames
            if args.numFrames > total_steps:
                # the user has asked for more frames than timesteps
                logging.warning("Number of frames --numFrames cannot be greater than number of timesteps --steps. Setting number of frames to number of steps.")
                args.numFrames = total_steps
                report_interval = int(total_steps/num_frames)
            elif args.numFrames == 0:
                # the user has asked us to run a simulation with 0 frames
                logging.warning("Number of frames --numFrames cannot be zero. Reverting to default of NUM_TIMESTEPS/REPORT_INTERVAL, where REPORT_INTERVAL takes the default value of 1000")
                report_interval = 1000
                num_frames = int(total_steps/report_interval)
            else:
                report_interval = int(total_steps/num_frames)
    else: # --reportInterval and --reportFrequency equal to each other but not equal to 1000
        if args.numFrames != -1:
            logging.warning("     --reportFrequency/--reportInterval take priority over --numFrames. Ignoring user-specified --numFrames.")
        report_interval = args.reportInterval
        num_frames = int(total_steps/report_interval)
    # check number of frames
    if num_frames == 0:
        # the user has asked us to run a simulation with 0 frames
        logging.warning("You have requested a simulation that contains less than 1 complete frame. Setting number of frames to min(1000, number of timesteps).")
        num_frames = min([1000,total_steps])
    assert num_frames <= total_steps, f"num_frames: {num_frames}, total_steps:{total_steps}" # this should be taken care of by earlier code
    # adjust total number of timesteps so that simulation ends on a complete frame
    if not (total_steps / num_frames).is_integer():
        logging.warning(f"    Number of timesteps --steps is not divisible by the number of frames. Increasing number of timesteps so that the simulation ends with a complete frame.")
        # we are not allowed to override the number of frames
        # we need to increase the value of total_steps such that it is an integer multiple of num_frames
        # (because we don't want to override the report frequency that we've already chosen)
        # we do this by rounding total_steps/num_frames up to the nearest integer, then multiplying by num_frames
        total_steps = (int(total_steps/num_frames)+1) * num_frames
        report_interval = total_steps / num_frames
    # make sure everything aggrees
    assert report_interval * num_frames == total_steps, f"report_interval: {report_interval}, num_frames: {num_frames}, total_steps: {total_steps}"
    assert report_interval.is_integer(), f"report_interval: {report_interval}, num_frames: {num_frames}, total_steps: {total_steps}"
    assert num_frames.is_integer(), f"report_interval: {report_interval}, num_frames: {num_frames}, total_steps: {total_steps}"
    report__interval = int(report_interval) # openmm functions want input type to be integer, not just integer-valued float
    num_frames = int(num_frames) # openmm functions want input type to be integer, not just integer-valued float
    total_steps = int(total_steps) # openmm functions want input type to be integer, not just integer-valued float
    # assign annealing parameters
    Tstart = args.tempStart
    Tend = args.tempEnd

    print(f"using force setup file from {forceSetupFile}")
    spec = importlib.util.spec_from_file_location("forces", forceSetupFile)
    # print(spec)
    forces = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(forces)


    oa = OpenMMAWSEMSystem(input_pdb_filename, k_awsem=1.0, chains=chain, xml_filename=openawsem.xml, seqFromPdb=seq, includeLigands=args.includeLigands)  # k_awsem is an overall scaling factor that will affect the relevant temperature scales
    myForces = forces.set_up_forces(oa, submode=args.subMode, contactParameterLocation=parametersLocation)
    # print(forces)
    # oa.addForces(myForces)

    if args.removeCMMotionRemover:
        oa.system.removeForce(0)
    oa.addForcesWithDefaultForceGroup(myForces)

    if args.fromCheckPoint:
        integrator = LangevinIntegrator(Tstart*kelvin, 1/picosecond, args.timeStep*femtoseconds)
        simulation = Simulation(oa.pdb.topology, oa.system, integrator, platform)
        simulation.loadCheckpoint(checkPointPath)
    else:
        # output the native and the structure after minimization
        integrator = CustomIntegrator(0.001)
        simulation = Simulation(oa.pdb.topology, oa.system, integrator, platform)
        simulation.context.setPositions(oa.pdb.positions)  # set the initial positions of the atoms
        simulation.reporters.append(PDBReporter(os.path.join(toPath, "native.pdb"), 1))
        simulation.reporters.append(DCDReporter(os.path.join(toPath, "movie.dcd"), 1))
        simulation.step(int(1))
        simulation.minimizeEnergy()  # first, minimize the energy to a local minimum to reduce any large forces that might be present
        simulation.step(int(1))


        # print("------------------Folding-------------------")
        # oa = OpenMMAWSEMSystem(input_pdb_filename, k_awsem=1.0, chains=chain, xml_filename=OPENAWSEM_LOCATION+"awsem.xml")  # k_awsem is an overall scaling factor that will affect the relevant temperature scales
        # myForces = forces.set_up_forces(oa, submode=args.subMode, contactParameterLocation=parametersLocation)
        # oa.addForces(myForces)

        integrator = LangevinIntegrator(Tstart*kelvin, 1/picosecond, args.timeStep*femtoseconds)
        # integrator.setRandomNumberSeed(A_NUMBER_AS_RANDOM_SEED)
        # integrator = CustomIntegrator(0.001)
        simulation = Simulation(oa.pdb.topology, oa.system, integrator, platform)
        # simulation.loadState(os.path.join(toPath, 'output.xml'))
        simulation.context.setPositions(oa.pdb.positions)  # set the initial positions of the atoms
        simulation.context.setVelocitiesToTemperature(Tstart*kelvin)  # set the initial velocities of the atoms according to the desired starting temperature
        # simulation.context.setVelocitiesToTemperature(Tstart*kelvin, A_RANDOM_SEED_NUMBER)
        simulation.minimizeEnergy()  # first, minimize the energy to a local minimum to reduce any large forces that might be present


    print("report_interval", report_interval)
    print("num_frames", num_frames)
    simulation.reporters.append(StateDataReporter(sys.stdout, report_interval, step=True, potentialEnergy=True, temperature=True))  # output energy and temperature during simulation
    simulation.reporters.append(StateDataReporter(os.path.join(toPath, "output.log"), report_interval, step=True, potentialEnergy=True, temperature=True)) # output energy and temperature to a file
    simulation.reporters.append(PDBReporter(os.path.join(toPath, "movie.pdb"), reportInterval=report_interval))  # output PDBs of simulated structures
    simulation.reporters.append(DCDReporter(os.path.join(toPath, "movie.dcd"), reportInterval=report_interval, append=True))  # output PDBs of simulated structures
    # simulation.reporters.append(DCDReporter(os.path.join(args.to, "movie.dcd"), 1))  # output PDBs of simulated structures
    # simulation.reporters.append(PDBReporter(os.path.join(args.to, "movie.pdb"), 1))  # output PDBs of simulated structures
    simulation.reporters.append(CheckpointReporter(os.path.join(toPath, checkpoint_file), checkpoint_reporter_frequency))  # save progress during the simulation

    if args.dryRun:
        if args.simulation_mode == 1: # test temperature setting
            deltaT = (Tend - Tstart) / num_frames
            for i in range(num_frames):
                integrator.setTemperature((Tstart + deltaT*i)*kelvin)
        raise SystemExit("Simulation configured successfully")

    print("Simulation Starts")
    start_time = time.time()

    if args.simulation_mode == 0:
        simulation.step(int(args.steps))
    elif args.simulation_mode == 1:
        deltaT = (Tend - Tstart) / num_frames
        for i in range(num_frames):
            integrator.setTemperature((Tstart + deltaT*i)*kelvin)
            simulation.step(report_interval) 

            # simulation.saveCheckpoint('step_%d.chk' % i)
            # simulation.context.setParameter("k_membrane", 0)
            # if i < snapShotCount/2:
            #     simulation.context.setParameter("k_membrane", (i % 2) * k_mem)
            #     simulation.context.setParameter("k_single_helix_orientation_bias", (i % 2) * k_single_helix_orientation_bias)
            # else:
            #     simulation.context.setParameter("k_membrane", k_mem)
            #     simulation.context.setParameter("k_single_helix_orientation_bias", k_single_helix_orientation_bias)

            # simulation.context.setParameter("k_membrane", (i)*(k_mem/snapShotCount))
            # simulation.context.setParameter("k_single_helix_orientation_bias", (i)*(k_single_helix_orientation_bias/snapShotCount))
            # print(simulation.context.getParameter("k_membrane"))


    # simulation.step(int(1e6))

    time_taken = time.time() - start_time  # time_taken is in seconds
    hours, rest = divmod(time_taken,3600)
    minutes, seconds = divmod(rest, 60)
    print(f"---{hours} hours {minutes} minutes {seconds} seconds ---")

    timeFile = os.path.join(toPath, "time.dat")
    with open(timeFile, "w") as out:
        out.write(str(time_taken)+"\n")

    # accompany with analysis run
    simulation = None
    time.sleep(10)
    os.chdir(pwd)
    print(os.getcwd())
    if args.fasta == "":
        analysis_fasta = ""
    else:
        analysis_fasta = f"--fasta {args.fasta}"
    if args.includeLigands:
        additional_cmd = "--includeLigands"
    else:
        additional_cmd = ""
    os.system(f"{sys.executable} mm_analyze.py {args.protein} -t {os.path.join(toPath, 'movie.dcd')} --subMode {args.subMode} -f {args.forces} {analysis_fasta} {additional_cmd} -c {chain}")


def main(args=None):
    parser = argparse.ArgumentParser(
        description="This is a python3 script to automatically copy the template file and run simulations")

    parser.add_argument("protein", help="The name of the protein")
    parser.add_argument("--name", default="simulation", help="Name of the simulation")
    parser.add_argument("--to", default="./", help="location of movie file")
    parser.add_argument("-c", "--chain", type=str, default="-1")
    parser.add_argument("-t", "--thread", type=int, default=-1, help="default is using all that is available")
    parser.add_argument("-p", "--platform", type=str, default="OpenCL", choices=["OpenCL", "CPU", "HIP", "Reference", "CUDA"], help="Platform to run the simulation.")
    parser.add_argument("-s", "--steps", type=float, default=1e7, help="step size, default 1e7")
    parser.add_argument("--tempStart", type=float, default=800, help="Starting temperature")
    parser.add_argument("--tempEnd", type=float, default=200, help="Ending temperature")
    parser.add_argument("--fromCheckPoint", type=str, default=None, help="The checkpoint file you want to start from")
    parser.add_argument("-m", "--simulation_mode", type=int, default=1,
                    help="default 1,\
                            0: constant temperature,\
                            1: temperature annealing")
    parser.add_argument("--subMode", type=int, default=-1)
    parser.add_argument("-f", "--forces", default="forces_setup.py")
    parser.add_argument("--parameters", default=None)
    parser.add_argument("-r", "--reportFrequency", type=int, default=1000, help="Frequency of timesteps measured against the reporter (record frame every N timesteps). Deprecated.")
    parser.add_argument("-I", "--reportInterval", type=int, default=1000, help="Frequency of timesteps measured against the reporter (record frame every N timesteps). Default value: 1000")
    parser.add_argument("--numFrames", type=int, default=-1, help="Number of frames to record. Timesteps will be 'wasted' at the end of simulation if number of timesteps is not divisible by number of frames. Default value: -1 (defer to --reportFrequency/--reportInterval)")
    parser.add_argument("--fromOpenMMPDB", action="store_true", default=False)
    parser.add_argument("--fasta", type=str, default="crystal_structure.fasta")
    parser.add_argument("--timeStep", type=int, default=2)
    parser.add_argument("--includeLigands", action="store_true", default=False)
    parser.add_argument('--device',default=0, help='OpenCL/CUDA device index')
    parser.add_argument('--removeCMMotionRemover', action="store_true", default=False, help='Removes CMMotionRemover. Recommended for periodic boundary conditions and membrane simulations')
    parser.add_argument('--dryRun',action="store_true",default=False,help="If True, quit before beginning simulation")
    
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    with open('commandline_args.txt', 'a') as f:
        f.write(' '.join(sys.argv))
        f.write('\n')
    print(' '.join(sys.argv))

    run(args)

if __name__=="__main__":
    main()
