import pytest
import logging

def numFrames_reportFrequency_part_of_mm_run(args):
    ## Configure number of timesteps and number of frames/recording interval
    # do some basic checks on user input
    if args.steps < 1:
        raise ValueError("Number of timesteps --steps must be a positive integer")
    else:
        total_steps = args.steps
    if args.numFrames > total_steps:
        # the user has asked for more frames than timesteps
        raise ValueError("Number of frames --numFrames cannot be greater than number of timesteps --steps")
    elif args.numFrames == 0:
        # the user has asked us to run a simulation with 0 frames
        raise ValueError("Number of frames --numFrames cannot be zero")
    if int(total_steps/args.reportFrequency) == 0: # int trucates so it works like the floor function for positive numbers
        # the user has asked us to run a simulation with 0 frames
        raise ValueError("Number of frames implied by --reportFrequency (calculated as NUM_TIMESTEPS/REPORT_FREQUENCY) must be at least 1")
    if args.reportFrequency != -1 and args.numFrames != -1:
        # The user has specified both --reportFrequency and --numFrames. We raise an error unless these directions are perfectly unambiguous.
        if not (total_steps/args.reportFrequency).is_integer():
            raise ValueError("Number of timesteps --steps is not divisible by the number of frames implied by --reportFrequency (calculated as NUM_TIMESTEPS/REPORT_FREQUENCY).")
        elif args.numFrames != total_steps/args.reportFrequency:
            raise ValueError("Number of frames from --numFrames and number of frames implied by --reportFrequency (calculated as NUM_TIMESTEPS/REPORT_FREQUENCY) disagree.\n\
                             Hint: you only need to give either --reportFrequency or --numFrames on the command line, not both.") 
    # assign number of frames, report frequency, and annealing parameters
    if args.reportFrequency == args.numFrames == -1:
        # the user has specified neither, so we revert to default
        logging.warning("    Not specified: (number of frames --numFrames OR frame reporting interval --reportFrequency), therefore reverting to default of 400 frames.")
        num_frames = 400
        if total_steps < 400:
            raise ValueError("Number of frames cannot be greater than number of timesteps --steps")
        if not (total_steps / num_frames).is_integer():
            logging.warning(f"    Number of timesteps --steps is not divisible by the number of frames. Increasing number of timesteps so that the simulation ends with a complete frame.")
            # we are not allowed to override the number of frames
            # so we need to increase the value of total_steps such that is it an integer multiple of 
            # (because we can't have a fractional reporting frequency)
            # we do this by round total_steps/num_frames up to the nearest integer, then multiplying by num_frames
            total_steps = (int(total_steps/num_frames)+1) * num_frames
        reporter_frequency = total_steps / num_frames
    elif args.numFrames != -1:
        # the user has specified --numFrames
        # if the user also specified --reportFrequency, then report frequency agrees with --numFrames (due to above check)
        num_frames = args.numFrames
        if not (total_steps / num_frames).is_integer():
            logging.warning(f"    Number of timesteps --steps is not divisible by the number of frames. Increasing number of timesteps so that the simulation ends with a complete frame.")
            # we are not allowed to override the number of frames
            # so we need to increase the value of total_steps such that it is an integer multiple of num_frames
            # (because we can't have a fractional reporting frequency)
            # we do this by rounding total_steps/num_frames up to the nearest integer, then multiplying by num_frames
            total_steps = (int(total_steps/num_frames)+1) * num_frames
        reporter_frequency = total_steps / num_frames
    elif args.reportFrequency != -1:
        # the user has specified --reportFrequency but not --numFrames
        reporter_frequency = args.reportFrequency
        if not (total_steps/reporter_frequency).is_integer():
            logging.warning(f"    Number of timesteps --steps is not divisible by the recording interval. Increasing number of timesteps so that the simulation ends with a complete frame.")
            # we are not allowed to override the reporting interval
            # so we need to increase the value of total_steps such that it is an integer multiple of reporter_frequency
            # (because we can't have a fractional number of frames)
            # we do this by rounding total_steps/reporter_frequency up to the nearest integer, then multiplying by reporter_frequency
            total_steps = (int(total_steps/reporter_frequency)+1) * reporter_frequency
        num_frames = total_steps / reporter_frequency
    else:
        raise AssertionError(f"Logical error in if-elif-elif, which should catch every case. frame: {args.numFrames}, frequency: {args.reportFrequency}")
    Tstart = args.tempStart
    Tend = args.tempEnd
    # make sure everything aggrees
    assert reporter_frequency * num_frames == total_steps, f"reporter_frequency: {reporter_frequency}, num_frames: {num_frames}, total_steps: {total_steps}"
    assert reporter_frequency.is_integer(), f"reporter_frequency: {reporter_frequency}, num_frames: {num_frames}, total_steps: {total_steps}"
    assert num_frames.is_integer(), f"reporter_frequency: {reporter_frequency}, num_frames: {num_frames}, total_steps: {total_steps}"
    reporter_frequency = int(reporter_frequency) # openmm functions want input type to be integer, not just integer-valued float
    num_frames = int(num_frames) # openmm functions want input type to be integer, not just integer-valued float
    total_steps = int(total_steps) # openmm functions want input type to be integer, not just integer-valued float

def test_numFrames_reportFrequency_handling():
    class Args:
        def __init__(self,arg_list):
            self.steps = arg_list[0]
            self.numFrames = arg_list[1]
            self.reportFrequency = arg_list[2]
            self.tempStart = 0 # don't care about this for test
            self.tempEnd = 0 # don't care about this for test

    args = Args([100,-1,-1]) # numFrames will default to 400, so we won't have enough timesteps
    with pytest.raises(ValueError):
        numFrames_reportFrequency_part_of_mm_run(args)

    args = Args([100,-1,99]) # this should work
    numFrames_reportFrequency_part_of_mm_run(args)

    args = Args([100,-1,100]) # this should work
    numFrames_reportFrequency_part_of_mm_run(args)

    args = Args([100,-1,101]) # this shouldn't work because we won't have a single timestep
    with pytest.raises(ValueError):
        numFrames_reportFrequency_part_of_mm_run(args)

    args = Args([100,99,-1]) # this should work
    numFrames_reportFrequency_part_of_mm_run(args)

    args = Args([100,100,-1]) # this should work
    numFrames_reportFrequency_part_of_mm_run(args)

    args = Args([100,101,-1]) # this shouldn't work because we won't have a single timestep
    with pytest.raises(ValueError):
        numFrames_reportFrequency_part_of_mm_run(args)

    args = Args([100,100,100]) # this shouldn't work because the given report interval and number of frames disagree
    with pytest.raises(ValueError):
        numFrames_reportFrequency_part_of_mm_run(args)

    args = Args([100,10,11]) # this shouldn't work because the given report interval and number of frames disagree
    with pytest.raises(ValueError):
        numFrames_reportFrequency_part_of_mm_run(args)

    args = Args([100,10,9]) # this shouldn't work because the given report interval and number of frames disagree
    with pytest.raises(ValueError):
        numFrames_reportFrequency_part_of_mm_run(args)

    args = Args([100,100,1]) # this should work
    numFrames_reportFrequency_part_of_mm_run(args)

    args = Args([100,1,100]) # this should work
    numFrames_reportFrequency_part_of_mm_run(args)

    args = Args([500,-1,-1]) # this should work
    numFrames_reportFrequency_part_of_mm_run(args)

    args = Args([500,101,-1]) # this should work
    numFrames_reportFrequency_part_of_mm_run(args)

    args = Args([500,-1,101]) # this should work
    numFrames_reportFrequency_part_of_mm_run(args)

if __name__ == "__main__":
    test_numFrames_reportFrequency_handling()